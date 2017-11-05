//! @file IonFlow.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/IonFlow.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"
#include "cantera/zdplaskin.h"

using namespace std;

namespace Cantera
{

IonFlow::IonFlow(IdealGasPhase* ph, size_t nsp, size_t points) :
    FreeFlame(ph, nsp, points),
    m_do_plasma(false),
    m_do_elec_heat(false),
    m_maxwellian_electron(false),
    m_stage(1),
    m_inletVoltage(0.0),
    m_outletVoltage(0.0),
    m_kElectron(npos),
    m_elec_num_density(1e17),
    m_elec_field(0.0),
    m_elec_frequency(0.0),
    m_plasma_multiplier(1.0),
    m_electron_multiplier(0.0)
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

    for (size_t i = 0; i < zdplaskinNSpecies(); i++) {
        char cstring[20];
        zdplaskinGetSpeciesName(cstring, &i);
        string speciesName(cstring);
        size_t k = m_thermo->speciesIndex(speciesName);
        m_plasmaSpeciesIndex.push_back(k);
        cout << speciesName << " " << k << endl;
    }

    // collision list
    vector<string> collision_list;
    collision_list.push_back("N2");
    collision_list.push_back("O2");
    collision_list.push_back("CH4");
    collision_list.push_back("CO2");
    collision_list.push_back("H2O");
    collision_list.push_back("H2");
    collision_list.push_back("CO");

    for (size_t i = 0; i < collision_list.size(); i++) {
        size_t k = m_thermo->speciesIndex(collision_list[i]);
        if (k != npos ) {
            m_kCollision.push_back(k);
            cout << collision_list[i] << " " << m_kCollision[i] << endl;
        }
    }

    // no bound for electric potential
    setBounds(c_offset_P, -1.0e20, 1.0e20);

    m_refiner->setActive(c_offset_P, false);
    m_mobility.resize(m_nsp*m_points);
    m_do_poisson.resize(m_points,false);
    m_do_velocity.resize(m_points,true);
    m_electronPower.resize(m_points, 0.0);
    m_electronTemperature.resize(m_points, 0.0);
    m_electronMobility.resize(m_points, 0.0);
    m_electronDiff.resize(m_points, 0.0);
    m_wdotPlasma.resize(zdplaskinNSpecies(),m_points,0.0);
}

void IonFlow::resize(size_t components, size_t points){
    StFlow::resize(components, points);
    m_mobility.resize(m_nsp*m_points);
    m_do_species.resize(m_nsp,true);
    m_do_poisson.resize(m_points,false);
    m_do_velocity.resize(m_points,true);
    m_fixedElecPoten.resize(m_points,0.0);
    m_fixedVelocity.resize(m_points);
    m_electronPower.resize(m_points, 0.0);
    m_electronTemperature.resize(m_points, 0.0);
    m_electronMobility.resize(m_points, 0.0);
    m_electronDiff.resize(m_points, 0.0);
    m_wdotPlasma.resize(zdplaskinNSpecies(),m_points,0.0);
}

void IonFlow::updateProperties(size_t jg, double* x, double* rsd, int* diag,
                              double rdt, size_t jmin, size_t jmax)
{
    StFlow::updateProperties(jg, x, rsd, diag, rdt, jmin, jmax);

    // update EEDF
    if (m_do_plasma) {
        for (size_t j = 0; j < m_points; j++) {
            for (size_t k : m_kCollision) {
                double number_density = ND(x,k,j);
                const char* species = m_thermo->speciesName(k).c_str();
                zdplaskinSetDensity(species, &number_density);
            }
            const double Tgas = T(x,j);
            double total_number_density = ND_t(j);
            zdplaskinSetGasTemp(&Tgas);
            zdplaskinSetElecField(&m_elec_field, &m_elec_frequency, &total_number_density);
            // get plasma properties
            m_electronTemperature[j] = zdplaskinGetElecTemp();
            m_electronMobility[j] = zdplaskinGetElecMobility(&total_number_density);
            m_electronDiff[j] = zdplaskinGetElecDiffCoeff();
            m_electronPower[j] = zdplaskinGetElecPower(&total_number_density);

            double* wdot_plasma = NULL;
            zdplaskinGetPlasmaSource(&wdot_plasma);
            for (size_t i = 0; i < zdplaskinNSpecies(); i++) {
                m_wdotPlasma(i,j) = wdot_plasma[i];
            }
        }
    } else {
        for (size_t j = 0; j < m_points; j++) {
            // get plasma properties
            m_electronTemperature[j] = T(x,j);
            m_electronMobility[j] = 0.4;
            m_electronDiff[j] = 0.4*(Boltzmann * T(x,j)) / ElectronCharge;
            m_electronPower[j] = 0.0;
        }
    }
}

void IonFlow::updateTransport(double* x, size_t j0, size_t j1)
{
    StFlow::updateTransport(x, j0, j1);
    for (size_t j = j0; j < j1; j++) {
        setGasAtMidpoint(x,j);
        m_trans->getMobilities(&m_mobility[j*m_nsp]);
        if (m_kElectron != npos) {
            size_t k = m_kElectron;
            m_mobility[k+m_nsp*j] = 0.5*(m_electronMobility[j]+m_electronMobility[j+1]);
            m_diff[k+m_nsp*j] = 0.5*(m_electronDiff[j]+m_electronDiff[j+1]);
        }
    }
}

void IonFlow::evalResidual(double* x, double* rsd, int* diag,
                           double rdt, size_t jmin, size_t jmax)
{
    StFlow::evalResidual(x, rsd, diag, rdt, jmin, jmax);
    if (m_stage == 3) {
        for (size_t j = jmin; j <= jmax; j++) {
            if (j == 0) {
                //rsd[index(c_offset_P, j)] = -phi(x,j);
                rsd[index(c_offset_P, j)] = E(x,j) - 0.0;
                diag[index(c_offset_P, j)] = 0;
            } else if (j == m_points - 1) {
                rsd[index(c_offset_P, j)] = m_outletVoltage - phi(x,j);
                //rsd[index(c_offset_P, j)] = E(x,j-1) - 0.0;
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

                // This method is used when you disable energy equation
                // but still maintain the velocity profile
                if (!m_do_velocity[j]) {
                    rsd[index(c_offset_U, j)] = u(x,j) - u_fixed(j);
                    diag[index(c_offset_U, j)] = 0;
                }
            }
        }
    }

    if (m_do_plasma) {
        for (size_t j = jmin; j <= jmax; j++) {
            if (j != 0 && j != m_points -1) {
                for (size_t i = 0; i < zdplaskinNSpecies(); i++) {
                    size_t k = m_plasmaSpeciesIndex[i];
                    if (k != npos) {
                        // multiply by the multiplier
                        m_wdotPlasma(i,j) *= m_plasma_multiplier;
                        rsd[index(c_offset_Y + k, j)] += m_wt[k] * m_wdotPlasma(i,j) / m_rho[j];
                    }
                }

                // update electron power
                if (m_do_elec_heat) {
                    if (m_do_energy[j]) {
                        rsd[index(c_offset_T, j)] += m_electronPower[j]
                                                     * abs(ND(x,m_kElectron,j))
                                                     / (m_rho[j] * m_cp[j])
                                                     * m_plasma_multiplier;
                    }
                }
            }
        }
    }
}

void IonFlow::updateDiffFluxes(const double* x, size_t j0, size_t j1)
{
    if (m_stage == 1) {
        frozenIonMethod(x,j0,j1);
    }
    if (m_stage == 2) {
        chargeNeutralityModel(x,j0,j1);
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

void IonFlow::chargeNeutralityModel(const double* x, size_t j0, size_t j1)
{
    for (size_t j = j0; j < j1; j++) {
        double wtm = m_wtm[j];
        double rho = density(j);
        double dz = z(j+1) - z(j);

        // mixture-average diffusion
        for (size_t k = 0; k < m_nsp; k++) {
            m_flux(k,j) = m_wt[k]*(rho*m_diff[k+m_nsp*j]/wtm);
            m_flux(k,j) *= (X(x,k,j) - X(x,k,j+1))/dz;
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

        // ambipolar diffusion
        double sum_chargeFlux = 0.0;
        double sum = 0.0;
        for (size_t k : m_kCharge) {
            double Xav = 0.5 * (X(x,k,j+1) + X(x,k,j));
            int q_k = m_speciesCharge[k];
            sum_chargeFlux += m_speciesCharge[k] / m_wt[k] * m_flux(k,j);
            // The mobility is used because it is more general than
            // using diffusion coefficient and Einstein relation
            sum += m_mobility[k+m_nsp*j] * Xav * q_k * q_k;
        }
        for (size_t k : m_kCharge) {
            double Xav = 0.5 * (X(x,k,j+1) + X(x,k,j));
            double drift;
            int q_k = m_speciesCharge[k];
            drift = q_k * q_k * m_mobility[k+m_nsp*j] * Xav / sum;
            drift *= -sum_chargeFlux * m_wt[k] / q_k;
            m_flux(k,j) += drift;
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

        // ambipolar diffusion
        double E_ambi = E(x,j);
        for (size_t k : m_kCharge) {
            double Yav = 0.5 * (Y(x,k,j) + Y(x,k,j+1));
            double drift = rho * Yav * E_ambi
                           * m_speciesCharge[k] * m_mobility[k+m_nsp*j];
            m_flux(k,j) += drift;
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

void IonFlow::setTransverseElecField(double elec_field, double elec_freq)
{
    m_elec_field = elec_field;
    m_elec_frequency = elec_freq;
}

void IonFlow::setPlasmaSourceMultiplier(double multiplier)
{
    m_plasma_multiplier = multiplier;
}

void IonFlow::setElectronTransportMultiplier(double multiplier)
{
    m_electron_multiplier = multiplier;
}

void IonFlow::enableElecHeat(bool withElecHeat)
{
    m_do_elec_heat = withElecHeat;
}

double IonFlow::getElecMobility(size_t j)
{
    return m_electronMobility[j];
}

double IonFlow::getElecDiffCoeff(size_t j)
{
    return m_electronDiff[j];
}

double IonFlow::getElecTemperature(size_t j)
{
    return m_electronTemperature[j];
}

double IonFlow::getElecCollisionHeat(size_t j)
{
    return m_electronPower[j];
}

void IonFlow::solvePlasma()
{
    bool changed = false;
    if (m_do_plasma) {
        changed = true;
    }
    m_do_plasma = true;
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    m_refiner->setActive(c_offset_P, true);
    if (changed) {
        needJacUpdate();
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

void IonFlow::solveVelocity(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_velocity[i]) {
                changed = true;
            }
            m_do_velocity[i] = true;
        }
    } else {
        if (!m_do_velocity[j]) {
            changed = true;
        }
        m_do_velocity[j] = true;
    }
    m_refiner->setActive(c_offset_U, true);
    m_refiner->setActive(c_offset_V, true);
    m_refiner->setActive(c_offset_T, true);
    if (changed) {
        needJacUpdate();
    }
}

void IonFlow::fixVelocity(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_velocity[i]) {
                changed = true;
            }
            m_do_velocity[i] = false;
        }
    } else {
        if (m_do_velocity[j]) {
            changed = true;
        }
        m_do_velocity[j] = false;
    }
    m_refiner->setActive(c_offset_U, false);
    m_refiner->setActive(c_offset_V, false);
    m_refiner->setActive(c_offset_T, false);
    if (changed) {
        needJacUpdate();
    }
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
    // save the velocity profile if the velocity is disabled
    bool v = m_do_velocity[0];
    for (size_t j = 0; j < m_points; j++) {
        if (!v) {
            m_fixedVelocity[j] = u(x,j);
        }
    }
    if (v) {
        solveVelocity();
    }
}

}
