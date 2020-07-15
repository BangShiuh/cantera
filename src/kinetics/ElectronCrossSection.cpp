/**
 *  @file ElectronCrossSection.cpp Definition file for class ElectronCrossSection.
 */
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ElectronCrossSection.h"

namespace Cantera {

ElectronCrossSection::ElectronCrossSection()
    : threshold(0.0)
{
}

ElectronCrossSection::~ElectronCrossSection()
{
}

void ElectronCrossSection::validate()
{
    if (kind == "ionization" || kind == "attachment" || kind == "excitation") {
        if (threshold < 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                               "The threshold of the process",
                               "(kind = '{}', target = '{}', product = '{}')",
                               "cannot be negative", kind, target, product);
        }
    } else if (kind != "effective" && kind != "elastic") {
        throw CanteraError("ElectronCrossSection::validate",
            "'{}' is an unknown type of cross section data.", kind);
    }
}

unique_ptr<ElectronCrossSection> newElectronCrossSection(const AnyMap& node)
{
    unique_ptr<ElectronCrossSection> ecs(new ElectronCrossSection());

    ecs->kind = node["kind"].asString();
    ecs->target = node["target"].asString();
    std::vector<vector_fp> data = node["data"].asVector<vector_fp>();
    // transpose data
    for (size_t i = 0; i < data.size(); i++) {
        ecs->energyLevel.push_back(data[i][0]);
        ecs->crossSection.push_back(data[i][1]);
    }
    ecs->threshold = node.getDouble("threshold", 0.0);
    ecs->product = node.getString("product", ecs->target);

    return ecs;
}

std::vector<shared_ptr<ElectronCrossSection>> getElectronCrossSections(const AnyValue& items)
{
    std::vector<shared_ptr<ElectronCrossSection> > all_cross_sections;
    for (const auto& node : items.asVector<AnyMap>()) {
        all_cross_sections.emplace_back(newElectronCrossSection(node));
    }
    return all_cross_sections;
}

}
