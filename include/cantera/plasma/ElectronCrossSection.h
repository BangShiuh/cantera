//! @file ElectronCrossSection.h Declaration for class Cantera::ElectronCrossSection.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ELECTRONCROSSSECTION_H
#define CT_ELECTRONCROSSSECTION_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

//! Contains data about the cross sections of electron collision
class ElectronCrossSection
{
public:
    ElectronCrossSection();

    //! ElectronCrossSection objects are not copyable or assignable
    ElectronCrossSection(const ElectronCrossSection&) = delete;
    ElectronCrossSection& operator=(const ElectronCrossSection& other) = delete;
    ~ElectronCrossSection();

    //! Validate the cross-section data.
    void validate();

    //! The name of the kind of electron collision
    std::string kind;

    //! The name of the target of electron collision
    std::string target;

    //! The product of electron collision
    std::string product;

    //! Data of cross section. data[i][j] where i is the index of a data point,
    //! j=0 is the electron energy [eV], and j=1 is the cross section [m^2].
    std::vector<vector_fp> data;

    //! The threshold of a process in [eV]
    double threshold;

    //! Extra data used for specific models
    AnyMap extra;
};

//! create an ElectronCrossSection object to store data.
unique_ptr<ElectronCrossSection> newElectronCrossSection(const AnyMap& node);

//! Get a vector of ElectronCrossSection objects to access the data.
std::vector<shared_ptr<ElectronCrossSection>> getElectronCrossSection(const AnyValue& items);

}

#endif
