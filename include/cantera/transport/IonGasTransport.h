/**
 * @file GasTransport.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_ION_GAS_TRANSPORT_H
#define CT_ION_GAS_TRANSPORT_H

#include "MixTransport.h"

namespace Cantera
{

class IonGasTransport : public MixTransport
{
public:
    //! Default constructor.
    IonGasTransport();

    virtual std::string transportType() const {
        return "Ion";
    }

protected:
    virtual void setupMM();

    // setup parameters for n64 model
    void setupN64();


    virtual void fitDiffCoeffs(MMCollisionInt& integrals);
    double omega11_n64(const double tstar, const double gamma);

    //! electrical properties
    vector_int m_speciesCharge;

    //! index of species with charges
    std::vector<size_t> m_kCharge;

    //! index of neutral species
    std::vector<size_t> m_kNeutral;

    //! index of electron
    size_t m_kElectron;

    DenseMatrix m_gamma;
};

}

#endif
