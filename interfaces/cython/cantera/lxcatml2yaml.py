#!/usr/bin/env python3
# encoding: utf-8

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
lxcatml2yaml.py: Convert the LXCat integral cross-section data in XML format (LXCATML) to YAML format

Usage:
    lxcatml2yaml [--input=<filename>]
                 [--database=<database name>]
                 [--species=<species name>]
                 [--output=<filename>]

Example:
    lxcatml2yaml --input=mycs.xml --database=itikawa --species=O2

If the output file name is not given, an output file with the same name as the
input file, with the extension changed to '.yaml'.
"""

from pathlib import Path
import argparse
import xml.etree.ElementTree as etree
from typing import Union
import sys
try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

BlockMap = yaml.comments.CommentedMap

# A class of yaml data for collision of a target species
class Process:
    def __init__(self, equation, energy_levels, cross_section):
        self.equation = equation
        self.energy_levels = energy_levels
        self.cross_section = cross_section

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('type', 'electron-collision-plasma'),
                        ('equation', node.equation),
                        ('energy-levels', node.energy_levels),
                        ('cross-section', node.cross_section),
                        ])
        return representer.represent_dict(out)

# Return indices of a child name
def get_children(parent, child_name):
    indices = []
    for i, child in enumerate(parent):
        if child.tag.find(child_name) != -1:
            indices.append(i)
    return [parent[x] for x in indices]

def FlowList(*args, **kwargs):
    """A YAML sequence that flows onto one line."""
    lst = yaml.comments.CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst

def FlowMap(*args, **kwargs):
    """A YAML mapping that flows onto one line."""
    m = yaml.comments.CommentedMap(*args, **kwargs)
    m.fa.set_flow_style()
    return m

def convert(
        inpfile: Union[str, Path] = None,
        databaseName: str = None,
        speciesName: str = None,
        outfile: Union[str, Path] = None,
        text: str = None
    ) -> None:
    """Convert an lxcat XML (LXCATML) file to a YAML file.

    :param inpfile:
        The input LXCATML file name. Exclusive with ``text``, only one of the two can be
        specified.
    :param outfile:
        The output YAML file name.
    :param databaseName:
        The name of the database. E.g. itikawa.
    :param speciesName:
        The target species name for electron collision. E.g. O2.
    :param text:
        Contains a string with the LXCATML input file content. Exclusive with ``inpfile``,
        only one of the two can be specified.

    All files are assumed to be relative to the current working directory of the Python
    process running this script.
    """
    if inpfile is not None and text is not None:
        raise ValueError("Only one of 'inpfile' or 'text' should be specified.")
    elif inpfile is not None:
        inpfile = Path(inpfile)
        lxcatml_text = inpfile.read_text().lstrip()
        if outfile is None:
            outfile = inpfile.with_suffix(".yaml")
    elif text is not None:
        if outfile is None:
            raise ValueError("If 'text' is passed, 'outfile' must also be passed.")
        lxcatml_text = text.lstrip()
    else:
        raise ValueError("One of 'inpfile' or 'text' must be specified")


    # Define yaml emitter
    emitter = yaml.YAML()
    emitter.register_class(Process)

    xml_tree = etree.fromstring(lxcatml_text)

    # Append all process together
    process_list = []

    for database in xml_tree:
        if database.attrib["id"] != databaseName:
            continue

        # Get groups node
        groups_node = get_children(database, "groups")[0]

        for group in groups_node:
            for process in get_children(group, "processes")[0]:
                # Target
                target = group.attrib["id"]
                if target != speciesName:
                    continue

                # Threshold
                threshold = 0.0
                parameters_node = get_children(process, "parameters")[0]
                if len(get_children(parameters_node, "parameter")) == 1:
                    parameter = get_children(parameters_node, "parameter")[0]
                    if parameter.attrib["name"] == 'E':
                        threshold = float(parameter.text)

                # Equation
                product_array=[]
                for product_node in get_children(process, "products")[0]:
                    if product_node.tag.find("electron") != -1:
                        product_array.append("e")
                    if product_node.tag.find("molecule") != -1:
                        product_name = product_node.text
                        if "state" in product_node.attrib:
                            state = product_node.attrib["state"]
                            product_name += f"({state})"
                        if "charge" in product_node.attrib:
                            charge = int(product_node.attrib["charge"])
                            if charge > 0:
                                product_name += charge*"+"
                            else:
                                product_name += -charge*"-"
                        product_array.append(product_name)

                products = " + ".join(product_array)
                equation = f"{target} + e => {products}"

                # Data
                data_x = get_children(process, "data_x")[0]
                data_y = get_children(process, "data_y")[0]

                energy_levels = FlowList(map(float, data_x.text.split(" ")))
                cross_section = FlowList(map(float, data_y.text.split(" ")))

                # edit energy levels and cross section
                if len(energy_levels) != len(cross_section):
                    raise ValueError("energy levels and cross section must have the same length")

                if energy_levels[0] > threshold:
                    energy_levels = FlowList([threshold, *energy_levels])
                    cross_section = FlowList([0.0, *cross_section])
                else:
                    cross_section[0] = 0.0

                # Save process
                process_list.append(Process(equation=equation,
                                            energy_levels=energy_levels,
                                            cross_section=cross_section))

    # Put process list in collision node
    collision_node = {"collisions": process_list}
    with Path(outfile).open("w") as output_file:
        emitter.dump(collision_node, output_file)

def main():
    """Parse command line arguments and pass them to `convert`."""
    parser = argparse.ArgumentParser(
        description="Convert the LXCat integral cross-section data in XML format (LXCATML) to YAML format",
        epilog=(
            "The 'output' argument is optional. If it is not given, an output "
            "file with the same name as the input file is used, with the extension "
            "changed to '.yaml'."
        ),
    )
    parser.add_argument("--input", required=True, type=str, help="The input LXCATML filename. Must be specified.")
    parser.add_argument("--database", required=True, type=str, help="The name of the database. Optional.")
    parser.add_argument("--species", required=True, type=str, help="The target species. Optional.")
    parser.add_argument("--output", nargs="?", help="The output YAML filename. Optional.")
    if len(sys.argv) not in [4, 5]:
        if len(sys.argv) > 5:
            print(
                "lxcatml2yaml.py: error: unrecognized arguments:",
                ' '.join(sys.argv[5:]),
                file=sys.stderr,
            )
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    input_file = Path(args.input)
    if args.output is None:
        output_file = input_file.with_suffix(".yaml")
    else:
        output_file = Path(args.output)

    convert(input_file, args.database, args.species, output_file)

if __name__ == "__main__":
    main()
