//! @file zdplaskin.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_ZDPLASKIN_H
#define CT_ZDPLASKIN_H

extern "C"
{
    void zdplaskinInit();
    void zdplaskinSetDensity(const char* cstring, const double* DENS);

    void zdplaskinSetConditions(const double* GAS_TEMPERATURE,
                                const double* REDUCED_FIELD);
    double zdplaskinGetElecTemp();
    double zdplaskinGetElecDiffCoeff();
    double zdplaskinGetElecMobility(const double*);
    void zdplaskinGetPlasmaSource(double** array);
    size_t zdplaskinNSpecies();
    void zdplaskinGetSpeciesName(char* cstring, size_t* index);
}

#endif

