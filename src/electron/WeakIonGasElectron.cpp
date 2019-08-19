// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#include "cantera/electron/WeakIonGasElectron.h"
#include "cantera/electron/ElectronFactory.h"
#include "cantera/base/utilities.h"
#include "cantera/numerics/funcs.h"
#include <Eigen/SparseLU>

namespace Cantera {

WeakIonGasElectron::WeakIonGasElectron()
    : m_chemionScatRate(0.0)
{
}

void WeakIonGasElectron::calculateTotalCrossSection()
{
    m_totalCrossSectionC.clear();
    m_totalCrossSectionC.resize(m_points, 0.0);
    m_totalCrossSectionB.clear();
    m_totalCrossSectionB.resize(m_points + 1, 0.0);
    for (size_t k = 0; k < m_ncs; k++) {
        vector_fp x = m_crossSections[k][0];
        vector_fp y = m_crossSections[k][1];
        for (size_t i = 0; i < m_points; i++) {
            m_totalCrossSectionC[i] += m_moleFractions[k] * linearInterp(m_gridC[i], x, y);
        }
        for (size_t i = 0; i < m_points + 1; i++) {
            m_totalCrossSectionB[i] += m_moleFractions[k] * linearInterp(m_gridB[i], x, y);
        }
    }
}

void WeakIonGasElectron::calculateTotalElasticCrossSection()
{
    m_sigmaElastic.clear();
    m_sigmaElastic.resize(m_points, 0.0);
    for (size_t k : m_kElastic) {
        vector_fp x = m_crossSections[k][0];
        vector_fp y = m_crossSections[k][1];
        for (size_t i = 0; i < m_points; i++) {
            m_sigmaElastic[i] += 2.0 * m_massRatios[k] * m_moleFractions[k] *
                                 linearInterp(m_gridB[i], x, y);
        }
    }
}

void WeakIonGasElectron::calculateDistributionFunction()
{
    // check if T or C is changed
    update_T();
    update_C();
    if (m_f0_ok == true) {
        return;
    }

    calculateTotalCrossSection();
    calculateTotalElasticCrossSection();

    double kT = m_kT;
    if (m_init_kTe != 0.0) {
        kT = m_init_kTe;
    }

    for (size_t j = 0; j < m_points; j++) {
        m_f0(j) = 2.0 * pow(1.0/Pi, 0.5) * pow(kT, -3./2.) *
                  std::exp(-m_gridC[j]/kT);
    }

    if (m_E != 0.0) {
        m_f0 = converge(m_f0);
        double Te = electronTemperature(m_f0);
        if (Te < m_thermo->temperature()) {
            for (size_t j = 0; j < m_points; j++) {
                m_f0(j) = 2.0 * pow(1.0/Pi, 0.5) * pow(m_kT, -3./2.) *
                          std::exp(-m_gridC[j]/m_kT);
            }
        }
    }
    m_f0_ok = true;
}

double WeakIonGasElectron::norm(Eigen::VectorXd f)
{
    vector_fp p(f.size());
    for (int i = 0; i < f.size(); i++) {
        p[i] = f(i) * pow(m_gridC[i], 0.5);
    }
    return simpsonQuadrature(m_gridC, p);
}

double WeakIonGasElectron::integralPQ(double a, double b, double u0, double u1,
                                      double g, double x0)
{
    double A1;
    double A2;
    if (g != 0.0) {
        double expm1a = std::expm1(g * (-a + x0));
        double expm1b = std::expm1(g * (-b + x0));
        double ag = a * g;
        double ag1 = ag + 1;
        double bg = b * g;
        double bg1 = bg + 1;
        A1 = (expm1a * ag1 + ag - expm1b * bg1 - bg) / (g*g);
        A2 = (expm1a * (2 * ag1 + ag * ag) + ag * (ag + 2) -
              expm1b * (2 * bg1 + bg * bg) - bg * (bg + 2)) / (g*g*g);
    } else {
        A1 = 0.5 * (b*b - a*a);
        A2 = 1.0 / 3.0 * (b*b*b - a*a*a);
    }

    // The interpolation formula of u(x) = c0 + c1 * x
    double c0 = (a * u1 - b * u0) / (a - b);
    double c1 = (u0 - u1) / (a - b);

    return c0 * A1 + c1 * A2;
}

vector_fp WeakIonGasElectron::vector_g(Eigen::VectorXd f0)
{
    vector_fp g(m_points, 0.0);
    g[0] = std::log(f0(1)/f0(0)) / (m_gridC[1] - m_gridC[0]);
    double N = m_points - 1;
    g[N] = std::log(f0(N)/f0(N-1)) / (m_gridC[N] - m_gridC[N-1]);
    for (size_t i = 1; i < m_points - 1; i++) {
        g[i] = std::log(f0(i+1)/f0(i-1)) / (m_gridC[i+1] - m_gridC[i-1]);
    }
    return g;
}

SpMat WeakIonGasElectron::matrix_PQ(vector_fp g, size_t k)
{
    std::vector<T> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        double i = m_i[k][n];
        double j = m_j[k][n];
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridC[j]);
        double q = m_inFactor[k] * m_gamma * m_moleFractions[k] * r;
        double p = - m_gamma * m_moleFractions[k] * r;
        tripletList.push_back(T(i, j, q));
        tripletList.push_back(T(j, j, p));
    }
    SpMat PQ(m_points, m_points);
    PQ.setFromTriplets(tripletList.begin(), tripletList.end());
    return PQ;
}

SpMat WeakIonGasElectron::matrix_P(vector_fp g, size_t k)
{
    std::vector<T> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        double j = m_j[k][n];
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridC[j]);
        double p = m_gamma * r;
        tripletList.push_back(T(j, j, p));
    }
    SpMat P(m_points, m_points);
    P.setFromTriplets(tripletList.begin(), tripletList.end());
    return P;
}

SpMat WeakIonGasElectron::matrix_Q(vector_fp g, size_t k)
{
    std::vector<T> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        double i = m_i[k][n];
        double j = m_j[k][n];
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridC[j]);
        double q = m_inFactor[k] * m_gamma * r;
        tripletList.push_back(T(i, j, q));
    }
    SpMat Q(m_points, m_points);
    Q.setFromTriplets(tripletList.begin(), tripletList.end());
    return Q;
}

SpMat WeakIonGasElectron::matrix_A(Eigen::VectorXd f0)
{
    vector_fp a0(m_points + 1);
    vector_fp a1(m_points + 1);
    size_t N = m_points - 1;
    // Scharfetter-Gummel scheme
    double nu = netProductionFreq(f0);
    a0[0] = NAN;
    a1[0] = NAN;
    a0[N+1] = NAN;
    a1[N+1] = NAN;
    for (size_t j = 1; j < m_points; j++) {
        double sigma_tilde = m_totalCrossSectionB[j] + nu / pow(m_gridB[j], 0.5) / m_gamma;
        double omega = 2 * Pi * m_F;
        double q = omega / (m_N * m_gamma * pow(m_gridB[j], 0.5));
        double W = -m_gamma * m_gridB[j] * m_gridB[j] * m_sigmaElastic[j];
        double F = sigma_tilde * sigma_tilde / (sigma_tilde * sigma_tilde + q * q);
        double DA = m_gamma / 3.0 * pow(m_E / m_N, 2.0) * m_gridB[j];
        double DB = m_gamma * m_kT * m_gridB[j] * m_gridB[j] * m_sigmaElastic[j];
        double D = DA / sigma_tilde * F + DB;
        double z = W * (m_gridC[j] - m_gridC[j-1]) / D;
        a0[j] = W / (1 - std::exp(-z));
        a1[j] = W / (1 - std::exp(z));
    }

    std::vector<T> tripletList;
    // center diagonal
    // zero flux b.c. at energy = 0
    tripletList.push_back(T(0, 0, a0[1]));
    for (size_t j = 1; j < m_points - 1; j++) {
        tripletList.push_back(T(j, j, a0[j+1] - a1[j]));
    }

    // upper diagonal
    for (size_t j = 0; j < m_points - 1; j++) {
        tripletList.push_back(T(j, j+1, a1[j+1]));
    }

    // lower diagonal
    for (size_t j = 1; j < m_points; j++) {
        tripletList.push_back(T(j, j-1, -a0[j]));
    }

    // zero flux b.c.
    tripletList.push_back(T(N, N, -a1[N]));

    SpMat A(m_points, m_points);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    //plus G
    SpMat G(m_points, m_points);
    for (size_t i = 0; i < m_points; i++) {
        G.insert(i,i) = 2.0 / 3.0 * (pow(m_gridB[i+1], 1.5) - pow(m_gridB[i], 1.5)) * nu;
    }
    return A + G;
}

Eigen::VectorXd WeakIonGasElectron::iterate(Eigen::VectorXd f0, double delta)
{
    SpMat PQ(m_points, m_points);
    vector_fp g = vector_g(f0);
    for (size_t k : m_kInelastic) {
        PQ += matrix_PQ(g, k);
    }

    SpMat A = matrix_A(f0);
    SpMat I(m_points, m_points);
    for (size_t i = 0; i < m_points; i++) {
        I.insert(i,i) = 1.0;
    }
    A -= PQ;
    A *= delta;
    A += I;

    // add chemionization scattering-in rate at the first grid
    f0(0) += m_chemionScatRate;

    // solve f0
    Eigen::SparseLU<SpMat> solver(A);
    Eigen::VectorXd f1 = solver.solve(f0);

    f1 *= 1.0 / norm(f1);
    return f1;
}

Eigen::VectorXd WeakIonGasElectron::converge(Eigen::VectorXd f0)
{
    double err0 = 0.0;
    double err1 = 0.0;
    double delta = m_delta0;
    for (size_t n = 0; n < m_maxn; n++) {
        if (0.0 < err1 && err1 < err0) {
            // log extrapolation attempting to reduce the error a factor m
            delta *= std::log(m_factorM) / (std::log(err0) - std::log(err1));
        }
        Eigen::VectorXd f1 = iterate(f0, delta);
        err0 = err1;
        Eigen::VectorXd Df0(m_points);
        for (size_t i = 0; i < m_points; i++) {
            Df0(i) = std::abs(f0(i) - f1(i));
        }
        err1 = norm(Df0);
        if (err1 < m_rtol) {
            return f1;
        }
        f0 = f1;
    }
    throw CanteraError("WeakIonGasElectron::converge", "Convergence failed");
}

double WeakIonGasElectron::netProductionFreq(Eigen::VectorXd f0)
{
    double nu = m_chemionScatRate;
    vector_fp g = vector_g(f0);

    for (size_t k = 0; k < m_ncs; k++) {
        if (m_kinds[k] == "IONIZATION" ||
            m_kinds[k] == "ATTACHMENT") {
            SpMat PQ = matrix_PQ(g, k);
            Eigen::VectorXd s = PQ * f0;
            for (size_t i = 0; i < m_points; i++) {
                nu += s[i];
            }
        }
    }
    return nu;
}

double WeakIonGasElectron::electronDiffusivity()
{
    calculateDistributionFunction();
    vector_fp y(m_points, 0.0);
    double nu = netProductionFreq(m_f0);
    for (size_t i = 0; i < m_points; i++) {
        if (m_gridC[i] != 0.0) {
            y[i] = m_gridC[i] * m_f0(i) /
                   (m_totalCrossSectionC[i] + nu / m_gamma / pow(m_gridC[i], 0.5));
        }
    }
    return 1./3. * m_gamma * simpsonQuadrature(m_gridC, y) / m_N;
}

double WeakIonGasElectron::electronMobility()
{
    calculateDistributionFunction();
    double nu = netProductionFreq(m_f0);
    vector_fp y(m_points + 1, 0.0);
    for (size_t i = 1; i < m_points; i++) {
        // calculate df0 at i-1/2
        double df0 = (m_f0(i) - m_f0(i-1)) / (m_gridC[i] - m_gridC[i-1]);
        if (m_gridB[i] != 0.0) {
            y[i] = m_gridB[i] * df0 /
                   (m_totalCrossSectionB[i] + nu / m_gamma / pow(m_gridB[i], 0.5));
        }
    }
    return -1./3. * m_gamma * simpsonQuadrature(m_gridB, y) / m_N;
}

double WeakIonGasElectron::realMobility()
{
    calculateDistributionFunction();
    double nu = netProductionFreq(m_f0);
    vector_fp y(m_points + 1, 0.0);
    for (size_t i = 1; i < m_points; i++) {
        // calculate df0 at i-1/2
        double df0 = (m_f0(i) - m_f0(i-1)) / (m_gridC[i] - m_gridC[i-1]);
        if (m_gridB[i] != 0.0) {
            double Q = m_totalCrossSectionB[i] + nu / m_gamma / pow(m_gridB[i], 0.5);
            double q = 2.0 * Pi * m_F / (m_N * m_gamma * pow(m_gridB[i], 0.5));
            y[i] = m_gridB[i] * Q / (Q * Q + q * q) * df0;
        }
    }
    return -1./3. * m_gamma * simpsonQuadrature(m_gridB, y) / m_N;
}

double WeakIonGasElectron::powerGain()
{
    if (m_F != 0.0) {
        return m_E * m_E * electronMobility();
    } else {
        return m_E * m_E * realMobility();
    }
}

double WeakIonGasElectron::elasticPowerLoss()
{
    calculateDistributionFunction();
    double sum = 0.0;
    vector_fp y(m_points + 1, 0.0);
    for (size_t i = 1; i < m_points; i++) {
        // calculate df0 at i-1/2
        double df0 = (m_f0(i) - m_f0(i-1)) / (m_gridC[i] - m_gridC[i-1]);
        double f0 = 0.5 * (m_f0(i-1) + m_f0(i));
        y[i] = m_sigmaElastic[i] * (m_gridB[i] * m_gridB[i] * f0 + m_kT *  df0);
    }
    sum += m_gamma * simpsonQuadrature(m_gridB, y);
    return sum * m_N;
}

double WeakIonGasElectron::totalCollisionFreq()
{
    double sum = 0.0;
    for (size_t k = 0; k < m_ncs; k++) {
        sum += m_moleFractions[k] * rateCoefficient(k);
    }
    return sum * m_N;
}

double WeakIonGasElectron::biMaxwellFraction(size_t k)
{
    double y0 = 1.0;
    double y1 = std::exp(-m_thresholds[k] / m_kT);
    y0 *= 1.0 / (y0 + y1);
    return y0;
}

double WeakIonGasElectron::inelasticPowerLoss()
{
    calculateDistributionFunction();
    double sum = 0.0;
    for (size_t k : m_kInelastic) {
        if (m_kinds[k] == "EXCITATION") {
            double y_low = biMaxwellFraction(k);
            double y_up = 1.0 - y_low;
            sum += m_thresholds[k] * m_moleFractions[k] *
                   (y_low * rateCoefficient(k) -
                    y_up * inverseRateCoefficient(k));
        }
    }
    return sum * m_N;
}

double WeakIonGasElectron::rateCoefficient(size_t k)
{
    calculateDistributionFunction();
    vector_fp g = vector_g(m_f0);
    SpMat P = matrix_P(g, k);
    Eigen::VectorXd s = P * m_f0;
    double sum = 0.0;
    for (size_t i = 0; i < m_points; i++) {
        sum += s(i);
    }
    return sum;
}

double WeakIonGasElectron::inverseRateCoefficient(size_t k)
{
    calculateDistributionFunction();
    vector_fp g = vector_g(m_f0);
    SpMat Q = matrix_Q(g, k);
    Eigen::VectorXd s = Q * m_f0;
    double sum = 0.0;
    for (size_t i = 0; i < m_points; i++) {
        sum += s(i);
    }
    return sum;
}

double WeakIonGasElectron::meanElectronEnergy()
{
    calculateDistributionFunction();
    double sum = 0;
    for (size_t i = 0; i < m_points - 1; i++) {
        sum += (pow(m_gridB[i+1], 2.5) -  pow(m_gridB[i], 2.5)) * m_f0(i);
    }
    return 0.4 * sum;
}

double WeakIonGasElectron::electronTemperature()
{
    double Te = 2./3. * meanElectronEnergy() / Boltzmann * ElectronCharge;
    if (Te < m_thermo->temperature()) {
        return m_thermo->temperature();
    } else {
        return Te;
    }
}

double WeakIonGasElectron::electronTemperature(Eigen::VectorXd f0)
{
    double sum = 0;
    for (size_t i = 0; i < m_points - 1; i++) {
        sum += (pow(m_gridB[i+1], 2.5) -  pow(m_gridB[i], 2.5)) * f0(i);
    }
    return 2./3. * 0.4 * sum / Boltzmann * ElectronCharge;
}

void WeakIonGasElectron::getNetPlasmaProductionRates(double* wdot)
{
    calculateDistributionFunction();
    for (size_t k = 0; k < m_thermo->nSpecies(); k++) {
        wdot[k] = 0.0;
    }
    size_t k_E = m_thermo->speciesIndex("E");
    if (k_E == npos) {
        return;
    }
    double X_E = m_thermo->moleFraction("E");
    for (size_t k : m_kInelastic) {
        size_t kr = m_thermo->speciesIndex(m_targets[k]);
        size_t kp0 = m_kProducts[k][0];
        size_t kp1 = m_kProducts[k][1];
        if (kr != npos && kp0 != npos && kp1 != npos) {
            double rate = m_N * m_N * m_moleFractions[k] * X_E *
                          rateCoefficient(k) / Avogadro;
            wdot[kp0] += rate;
            if (kp1 != 9999) {
                wdot[kp1] += rate;
            }
            wdot[kr] -= rate;
            wdot[k_E] += (m_inFactor[k] - 1) * rate;
        }
    }
}

}
