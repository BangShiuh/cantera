# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

import warnings
import weakref

cdef class Collision:
    """
 
    """
    def __cinit__(self, *args, init=True, **kwargs):
        if init:
            self._collision.reset(new CxxCollision())
            self.collision = self._collision.get()

    def __init__(self, kind=None, threshold=None, *args, init=True, **kwargs):
        if not init:
            return

        if threshold is not None:
            self.collision.threshold = threshold

        if kind is not None:
            self.collision.kind = stringify(kind)

    cdef _assign(self, shared_ptr[CxxCollision] other):
        self._collision = other
        self.collision = self._collision.get()

    @staticmethod
    def fromYaml(text, Kinetics kinetics):
        """
        Create a Collision object from its YAML string representation.
        """
        cxx_collision = CxxNewCollision(AnyMapFromYamlString(stringify(text)),
                                        deref(kinetics.kinetics))
        collision = Collision(init=False)
        collision._assign(cxx_collision)
        return collision

    @staticmethod
    def listFromYaml(text, Kinetics kinetics, section=None):
        """
        Create a list of Collision objects from all the collisions defined in a YAML
        string. If ``text`` is a YAML mapping, the ``section`` name of the list
        to be read must be specified. If ``text`` is a YAML list, no ``section``
        name should be supplied.
        """
        root = AnyMapFromYamlString(stringify(text))

        # ``items`` is the pseudo-key used to access a list when it is at the
        # top level of a YAML document
        cxx_collisions = CxxGetCollisions(root[stringify(section or "items")],
                                          deref(kinetics.kinetics))
        collisions = []
        for a in cxx_collisions:
            b = Collision(init=False)
            b._assign(a)
            collisions.append(b)
        return collisions

    @staticmethod
    def listFromFile(filename, Kinetics kinetics, section='collisions'):
        """
        Create a list of Species objects from all of the species defined in a
        YAML file. Return species from the section *section*.

        Directories on Cantera's input file path will be searched for the
        specified file.
        """
        root = AnyMapFromYamlFile(stringify(filename))
        cxx_collisions = CxxGetCollisions(root[stringify(section)], deref(kinetics.kinetics))

        collisions = []
        for a in cxx_collisions:
            b = Collision(init=False)
            b._assign(a)
            collisions.append(b)
        return collisions

    property equation:
        """
        A string giving the chemical equation for this collision. Determined
        automatically based on `reactants` and `products`.
        """
        def __get__(self):
            return pystr(self.collision.equation)

    property kind:
        """
        The kind of the collision
        """
        def __get__(self):
            return pystr(self.collision.kind)

    property threshold:
        """
        The threshold of the collision in eV.
        """
        def __get__(self):
            return self.collision.threshold

    property energy_data:
        """
        The energy data used for interpolate the cross section
        """
        def __get__(self):
            return self.input_data["energy-data"]

    property cross_section_data:
        """
        The cross section data
        """
        def __get__(self):
            return self.input_data["cross-section-data"]

    property input_data:
        def __get__(self):
            cdef CxxThermoPhase* phase = self._phase.thermo if self._phase else NULL
            return anymap_to_dict(self.collision.parameters(phase))

    def update_user_data(self, data):
        """
        Add the contents of the provided `dict` as additional fields when generating
        YAML phase definition files with `Solution.write_yaml` or in the data returned
        by `input_data`. Existing keys with matching names are overwritten.
        """
        self.collision.input.update(dict_to_anymap(data), False)

    def clear_user_data(self):
        """
        Clear all saved input data, so that the data given by `input_data` or
        `Solution.write_yaml` will only include values generated by Cantera based on
        the current object state.
        """
        self.collision.input.clear()

    def __repr__(self):
        return f"<{self.__class__.__name__}: {self.equation}>"
