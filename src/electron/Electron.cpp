// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/electron/Electron.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include "cantera/numerics/funcs.h"
#include <iostream>
#include <limits>
#include <set>

namespace Cantera {

Electron::Electron()
    : m_electronCrossSectionTargets(0)
    , m_electronCrossSectionKinds(0)
    , m_ncs(0)
    , m_points(1000)
    , m_kTe(Undef)
    , m_kT(Undef)
    , m_electronCrossSections_ok(false)
    , m_f0_ok(false)
{
    // default energy grid
    m_eps.resize(m_points);
    for (size_t j = 0; j < m_points; j++) {
        m_eps[j] = j / 100.0;
    }

    m_gamma = pow(2.0 * ElectronCharge / ElectronMass, 0.5);
}

Electron::~Electron()
{
}

void Electron::init(thermo_t* thermo)
{
    m_thermo = thermo;
    m_f0_ok = false;
}

void Electron::update_T()
{
    // signal that temperature-dependent quantities will need to be recomputed
    // before use, and update the local temperature.
    m_kT = Boltzmann * m_thermo->temperature() / ElectronCharge;

    // flag for quantities need to be re-calculated
    m_f0_ok = false;
}

void Electron::update_C()
{
    // signal that concentration-dependent quantities will need to be recomputed
    // before use, and update the local mole fractions.
    compositionMap gas_composition = m_thermo->getMoleFractionsByName(0.0);
    m_moleFractions.resize(m_ncs, 0.0);
    for (auto const& x : gas_composition) {
        bool not_found = true;
        for (size_t i = 0; i < m_ncs; i++) {
            if (m_electronCrossSectionTargets[i] == x.first) {
                m_moleFractions[i] = x.second;
                not_found = false;
            }
        }
        if (not_found) {
            if (x.second > 0.01) {
                std::cout << "The mole fraction of species " << x.first
                            << " is more than 0.01 but it has no cross section data."
                            << std::endl;
            }
        }
    }
    // flag for quantities need to be re-calculated
    m_f0_ok = false;
}

bool Electron::addElectronCrossSection(shared_ptr<ElectronCrossSection> ecs)
{
    if (std::find(m_electronCrossSectionTargets.begin(),
                  m_electronCrossSectionTargets.end(),
                  ecs->target) != m_electronCrossSectionTargets.end()) {
        if (std::find(m_electronCrossSectionKinds.begin(),
                      m_electronCrossSectionKinds.end(),
                      ecs->kind) != m_electronCrossSectionKinds.end()) {
            throw CanteraError("Electron::addElectronCrossSection",
                                "Already contains a data of type '{}' for '{}'.",
                                ecs->kind, ecs->target);
        }
    }
    ecs->validate();
    m_electronCrossSectionTargets.push_back(ecs->target);
    m_electronCrossSectionKinds.push_back(ecs->kind);
    m_massRatios.push_back(ecs->mass_ratio);

    // transpose data
    std::vector<std::vector<double>> transdata(2, std::vector<double>(ecs->data.size()));
    for (size_t i = 0; i < ecs->data.size(); i++) {
        for (size_t j = 0; j < 2; j++) {
            transdata[j][i] = ecs->data[i][j];
        }
    }
    m_electronCrossSectionData.push_back(transdata);
    m_ncs++;
    m_electronCrossSections_ok = false;
    return true;
}

void Electron::setupGrid(size_t n, const double* eps)
{
    m_points = n;
    m_eps.resize(n);
    for (size_t j = 0; j < m_points; j++) {
        m_eps[j] = eps[j];
    }
    m_electronCrossSections_ok = false;
    m_f0_ok = false;
}

void Electron::setupCrossSections()
{
    m_electronCrossSections.resize(m_ncs, std::vector<double>(m_points));
    for (size_t i = 0; i < m_ncs; i++) {
        for (size_t j = 0; j < m_points; j++) {
            m_electronCrossSections[i][j] = linearInterp(m_eps[j],
                                                        m_electronCrossSectionData[i][0],
                                                        m_electronCrossSectionData[i][1]);
        }
    }
    m_electronCrossSections_ok = true;
}

void Electron::calculateTotalCrossSection()
{
    if (m_electronCrossSections_ok == false) {
        setupCrossSections();
    }
    m_totalCrossSection.resize(m_points);
    for (size_t i = 0; i < m_ncs; i++) {
        if (m_electronCrossSectionKinds[i] == "EFFECTIVE") {
            for (size_t j = 0; j < m_points; j++) {
                m_totalCrossSection[j] += m_moleFractions[i] * m_electronCrossSections[i][j];
            }
        }
    }
}

void Electron::calculateDistributionFunction()
{
    update_T();
    update_C();
    calculateTotalCrossSection();
    m_f0.resize(m_points);
    m_df0.resize(m_points);
    m_f0_ok = true;
}

double Electron::electronDiffusivity(double N)
{
    if (m_f0_ok == false) {
        calculateDistributionFunction();
    }
    double sum = 0.0;
    for (size_t j = 0; j < m_points - 1; j++) {
        sum += 0.5 * (m_eps[j] * m_f0[j] / m_totalCrossSection[j] +
                m_eps[j+1] * m_f0[j+1] / m_totalCrossSection[j+1]) *
                (m_eps[j+1] - m_eps[j]);
    }
    return 1./3. * m_gamma * sum / N;
}

double Electron::electronMobility(double N)
{
    if (m_f0_ok == false) {
        calculateDistributionFunction();
    }
    double sum = 0.0;
    for (size_t j = 0; j < m_points - 1; j++) {
        sum += 0.5 * (m_eps[j] * m_df0[j] / m_totalCrossSection[j] +
                m_eps[j+1] * m_df0[j+1] / m_totalCrossSection[j+1]) *
                (m_eps[j+1] - m_eps[j]);
    }
    return -1./3. * m_gamma * sum / N;
}

}
