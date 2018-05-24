//! @file IonFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/IonFlow.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"
#include "cantera/numerics/polyfit.h"

using namespace std;

namespace Cantera
{

IonFlow::IonFlow(IdealGasPhase* ph, size_t nsp, size_t points) :
    FreeFlame(ph, nsp, points),
    m_import_electron_transport(false),
    m_ohmic_heat_E(0.0),
    m_plasma_multiplier(0.0),
    m_stage(1),
    m_inletVoltage(0.0),
    m_outletVoltage(0.0),
    m_kElectron(npos)
{
    // make a local copy of species charge
    for (size_t k = 0; k < m_nsp; k++) {
        m_speciesCharge.push_back(m_thermo->charge(k));
    }

    // Find indices for charge of species
    for (size_t k = 0; k < m_nsp; k++){
        if (m_speciesCharge[k] != 0){
            m_kCharge.push_back(k);
        } else {
            m_kNeutral.push_back(k);
        }
    }

    // Find the index of electron
    if (m_thermo->speciesIndex("E") != npos ) {
        m_kElectron = m_thermo->speciesIndex("E");
    }

    // no bound for electric potential
    setBounds(c_offset_P, -1.0e20, 1.0e20);

    m_refiner->setActive(c_offset_P, false);
    m_mobility.resize(m_nsp*m_points);
    m_do_poisson.resize(m_points,false);
}

void IonFlow::resize(size_t components, size_t points){
    StFlow::resize(components, points);
    m_mobility.resize(m_nsp*m_points);
    m_do_species.resize(m_nsp,true);
    m_do_poisson.resize(m_points,false);
    m_fixedElecPoten.resize(m_points,0.0);
    m_elecMobility.resize(m_points);
    m_elecDiffCoeff.resize(m_points);
}

void IonFlow::updateTransport(double* x, size_t j0, size_t j1)
{
    StFlow::updateTransport(x,j0,j1);
    for (size_t j = j0; j < j1; j++) {
        setGasAtMidpoint(x,j);
        m_trans->getMobilities(&m_mobility[j*m_nsp]);
        double tlog = log(m_thermo->temperature());
        if (m_import_electron_transport) {
            size_t k = m_kElectron;
            m_mobility[k+m_nsp*j] *= (1.0 - m_plasma_multiplier);
            m_mobility[k+m_nsp*j] += m_plasma_multiplier * poly5(tlog, m_mobi_e_fix.data());
            m_diff[k+m_nsp*j] *= (1.0 - m_plasma_multiplier);
            m_diff[k+m_nsp*j] += m_plasma_multiplier * poly5(tlog, m_diff_e_fix.data());
        }
    }
}

void IonFlow::updateDiffFluxes(const double* x, size_t j0, size_t j1)
{
    if (m_stage == 1) {
        frozenIonMethod(x,j0,j1);
    }
    if (m_stage == 3) {
        poissonEqnMethod(x,j0,j1);
    }
}

void IonFlow::frozenIonMethod(const double* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double wtm = m_wtm[j];
        double rho = density(j);
        double dz = z(j+1) - z(j);
        double sum = 0.0;
        for (size_t k : m_kNeutral) {
            m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
            sum -= m_flux(k,j);
        }

        // correction flux to insure that \sum_k Y_k V_k = 0.
        for (size_t k : m_kNeutral) {
            m_flux(k,j) += sum*Y(x,k,j);
        }

        // flux for ions
        // Set flux to zero to prevent some fast charged species (e.g. electron)
        // to run away
        for (size_t k : m_kCharge) {
            m_flux(k,j) = 0;
        }
    }
}

void IonFlow::poissonEqnMethod(const double* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double wtm = m_wtm[j];
        double rho = density(j);
        double dz = z(j+1) - z(j);

        // mixture-average diffusion
        double sum = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
            sum -= m_flux(k,j);
        }

        // ambipolar diffusion
        double E_ambi = E(x,j);
        for (size_t k : m_kCharge) {
            double Yav = 0.5 * (Y(x,k,j) + Y(x,k,j+1));
            double drift = rho * Yav * E_ambi
                           * m_speciesCharge[k] * m_mobility[k+m_nsp*j];
            m_flux(k,j) += drift;
        }

        // correction flux
        double sum_flux = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum_flux -= m_flux(k,j); // total net flux
        }
        double sum_ion = 0.0;
        for (size_t k : m_kCharge) {
            sum_ion += Y(x,k,j);
        }
        // The portion of correction for ions is taken off
        for (size_t k : m_kNeutral) {
            m_flux(k,j) += Y(x,k,j) / (1-sum_ion) * sum_flux;
        }
    }
}

void IonFlow::setSolvingStage(const size_t stage)
{
    if (stage == 1 || stage == 2 || stage == 3) {
        m_stage = stage;
    } else {
        throw CanteraError("IonFlow::updateDiffFluxes",
                    "solution phase must be set to:"
                    "1: frozenIonMethod"
                    "2: chargeNeutralityModel"
                    "3: poissonEqnMethod");
    }
}

void IonFlow::setElectricPotential(const double v1, const double v2)
{
    // This method can be used when you want to add external voltage
    m_inletVoltage = v1;
    m_outletVoltage = v2;
}

void IonFlow::getWdot(double* x, size_t j)
{
    StFlow::getWdot(x,j);
    if (m_plasmaRateCoeff.size() > 0) {
        // index
        size_t k1 = m_thermo->speciesIndex("O2");
        size_t k2 = m_thermo->speciesIndex("O2(a1dg)");
        // obtain rate
        double tlog = log(T(x,j));
        double rate = poly5(tlog, m_plasmaRateCoeff.data());
        rate *= ND(x,k1,j) * ND(x,m_kElectron,j);
        rate /= Avogadro;
        m_wdot(k1,j) -= rate;
        m_wdot(k2,j) += rate;
    }
    if (m_electronTemperature.size() > 0) {
        double Tg = T(x,j);
        double Te = poly5(Tg, m_electronTemperature.data());
    }
}


void IonFlow::evalResidual(double* x, double* rsd, int* diag,
                           double rdt, size_t jmin, size_t jmax)
{
    StFlow::evalResidual(x, rsd, diag, rdt, jmin, jmax);
    for (size_t j = jmin; j <= jmax; j++) {
        if (j != 0 && j != m_points - 1) {
            // update electron power
            if (m_ohmic_heat_E > 0.0) {
                size_t k = m_kElectron;
                double mu_av = 0.5 * (m_mobility[k+m_nsp*j] + m_mobility[k+m_nsp*(j-1)]);
                if (m_do_energy[j]) {
                    double ND_e = (ND(x,k,j) > 0.0 ? ND(x,k,j) : 0.0);
                    rsd[index(c_offset_T, j)] += ElectronCharge * mu_av
                                                 * m_ohmic_heat_E * m_ohmic_heat_E
                                                 * ND_e / (m_rho[j] * m_cp[j]);
                }
            }
        }
    }
    if (m_stage == 3) {
        for (size_t j = jmin; j <= jmax; j++) {
            if (j == 0) {
                // enforcing the flux for charged species is difficult
                // since charged species are also affected by electric
                // force, so Neumann boundary condition is used.
                for (size_t k : m_kCharge) {
                    rsd[index(c_offset_Y + k, 0)] = Y(x,k,0) - Y(x,k,1);
                }
                rsd[index(c_offset_P, j)] = m_inletVoltage - phi(x,j);
                diag[index(c_offset_P, j)] = 0;
            } else if (j == m_points - 1) {
                rsd[index(c_offset_P, j)] = m_outletVoltage - phi(x,j);
                diag[index(c_offset_P, j)] = 0;
            } else {
                //-----------------------------------------------
                //    Poisson's equation
                //
                //    dE/dz = e/eps_0 * sum(q_k*n_k)
                //
                //    E = -dV/dz
                //-----------------------------------------------
                double chargeDensity = 0.0;
                for (size_t k : m_kCharge) {
                    chargeDensity += m_speciesCharge[k] * ElectronCharge * ND(x,k,j);
                }
                rsd[index(c_offset_P, j)] = dEdz(x,j) - chargeDensity / epsilon_0;
                diag[index(c_offset_P, j)] = 0;
            }
        }
    }
}

void IonFlow::solvePoissonEqn(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_poisson[i]) {
                changed = true;
            }
            m_do_poisson[i] = true;
        }
    } else {
        if (!m_do_poisson[j]) {
            changed = true;
        }
        m_do_poisson[j] = true;
    }
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    m_refiner->setActive(c_offset_P, true);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::fixElectricPotential(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_poisson[i]) {
                changed = true;
            }
            m_do_poisson[i] = false;
        }
    } else {
        if (m_do_poisson[j]) {
            changed = true;
        }
        m_do_poisson[j] = false;
    }
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    m_refiner->setActive(c_offset_P, false);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::setElectronTransport(vector_fp& tfix, vector_fp& diff_e,
                                   vector_fp& mobi_e)
{
    m_import_electron_transport = true;
    size_t degree = 5;
    size_t n = tfix.size();
    vector_fp tlog;
    for (size_t i = 0; i < n; i++) {
        tlog.push_back(log(tfix[i]));
    }
    vector_fp w(n, -1.0);
    m_diff_e_fix.resize(degree + 1);
    m_mobi_e_fix.resize(degree + 1);
    polyfit(n, degree, tlog.data(), diff_e.data(), w.data(), m_diff_e_fix.data());
    polyfit(n, degree, tlog.data(), mobi_e.data(), w.data(), m_mobi_e_fix.data());
}

void IonFlow::setPlasmaRateCoeff(vector_fp& tfix, vector_fp& k)
{
    size_t degree = 5;
    size_t n = tfix.size();
    vector_fp tlog;
    for (size_t i = 0; i < n; i++) {
        tlog.push_back(log(tfix[i]));
    }
    vector_fp w(n, -1.0);
    m_plasmaRateCoeff.resize(degree + 1);
    polyfit(n, degree, tlog.data(), k.data(), w.data(), m_plasmaRateCoeff.data());
}

void IonFlow::setElectronTemperature(vector_fp& tfix, vector_fp& Te)
{
    size_t degree = 5;
    size_t n = tfix.size();
    vector_fp w(n, -1.0);
    m_electronTemperature.resize(degree + 1);
    polyfit(n, degree, tfix.data(), Te.data(), w.data(), m_electronTemperature.data());
}

void IonFlow::setOhmicHeatingElectricField(const double efield)
{
    m_ohmic_heat_E = efield;
}

void IonFlow::setPlasmaMultiplier(const double multi)
{
    m_plasma_multiplier = multi;
}

void IonFlow::_finalize(const double* x)
{
    FreeFlame::_finalize(x);

    bool p = m_do_poisson[0];
    for (size_t j = 0; j < m_points; j++) {
        if (!p) {
            m_fixedElecPoten[j] = phi(x, j);
        }
    }
    if (p) {
        solvePoissonEqn();
    }
}

}
