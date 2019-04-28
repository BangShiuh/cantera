/**
 * @file WeakIonGasElectron.h
 * Header file for class WeakIonGasElectron.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_WEAKIONGASELECTRON_H
#define CT_WEAKIONGASELECTRON_H

#include "cantera/electron/Electron.h"
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

namespace Cantera
{
/**
 * This class calculates the properties of electron in a weakly ionized gas.
 * Only electron-neutral collisions are considered for calculating the
 * electron energy distribution function (EEDF).
 *
 * Reference:
 * [1] G. J. M. Hagelaar and L. C. Pitchford
 * "Solving the Boltzmann equation to obtain electron transport
 * coefficients and rate coefficients for fluid models."
 * Plasma Sources Science and Technology 14.4 (2005): 722.
 * doi: https://doi.org/10.1088/0963-0252/14/4/011
 * [2] A. Luque, "BOLOS: An open source solver for the Boltzmann equation,"
 * https://github.com/aluque/bolos.
 * [3] D. McElroy, C. Walsh, A. Markwick, M. Cordiner, K. Smith, T. Millar,
 * "The umist database for astrochemistry 2012,"
 * Astronomy & Astrophysics 550 (2013) A36.
 * doi: https://doi.org/10.1051/0004-6361/201220465
 * @ingroup electron
 */
class WeakIonGasElectron: public Electron
{
public:
    WeakIonGasElectron();

    virtual double electronDiffusivity(double N);
    virtual double electronMobility(double N);
    virtual double meanElectronEnergy();

    //! electron temperature
    //! If the reduced electric field is set, electron tempeature is calculated
    //! from EEDF.
    virtual double electronTemperature();

    //! Set chemionization scattering-in rate
    virtual void setChemionScatRate(double rate) {
        m_chemionScatRate = rate;
        m_setChemionScatRate = true;
        m_f0_ok = false;
    }

protected:
    //! Calculate distribution function
    void calculateDistributionFunction();

    //! Calculate total cross section
    void calculateTotalCrossSection();

    //! Calculate total elastic cross section
    void calculateTotalElasticCrossSection();

    //! Calculate chemionization scattering-in rate
    //! The reaction rate coefficient is from [3], and can be downloaded at
    //! http://udfa.ajmarkwick.net/
    void calculateChemionScatRate();

    //! Set grid cache
    void setGridCache();

    //! The integral in [a, b] of x * u(x) exp[g * (x0 - x)]
    //! assuming that u is linear with u(a) = u0 and u(b) = u1
    double integralPQ(double a, double b, double u0, double u1,
                       double g, double x0);

    //! Norm of electron energy distribution function
    double norm(Eigen::VectorXd f);

    //! Vector g is used by matrix_PQ.
    vector_fp vector_g(Eigen::VectorXd f0);

    //! The matrix of scattering-out and scattering-in
    SpMat matrix_PQ(Eigen::VectorXd f0, vector_fp g, size_t k);

    //! matrix A represents equation (45) of ref. [1]
    SpMat matrix_A(Eigen::VectorXd f0);

    //! An iteration of solving electron energy distribution function
    Eigen::VectorXd iterate(Eigen::VectorXd f0, double delta = 1e14);

    //! Iterate until convergence and obtain EEDF
    Eigen::VectorXd converge(Eigen::VectorXd f0);

    //! Reduced net production frequency. Equation (10) of ref. [1]
    //! divided by N.
    //! @param f0 EEDF
    //! @param chem_rate Chemionization rate. This is an optional source of
    //!        electron.
    double netProductionFreq(Eigen::VectorXd f0);

    //! Total electron cross section on the cell center of energy grid
    vector_fp m_totalCrossSectionC;

    //! Total electron cross section on the cell boundary (i-1/2) of
    //! energy grid
    vector_fp m_totalCrossSectionB;

    //! vector of total elastic cross section weighted with mass ratio
    vector_fp m_sigmaElastic;

    //! Location of cell j for grid cache
    std::vector<std::vector<size_t>> m_j;

    //! Location of cell i for grid cache
    std::vector<std::vector<size_t>> m_i;

    //! Cross section at the boundaries of the overlap of cell i and j
    std::vector<std::vector<vector_fp>> m_sigma;

    //! The energy boundaries of the overlap of cell i and j
    std::vector<std::vector<vector_fp>> m_eps;

    //! Set chemionization scattering-in rate
    double m_chemionScatRate;

    //! flag of setChemionScatRate
    bool m_setChemionScatRate;
};

}

#endif
