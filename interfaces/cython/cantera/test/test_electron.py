import numpy as np

import cantera as ct
from . import utilities
import copy


class TestElectron(utilities.CanteraTest):
    def setUp(self):
        self.gas = ct.Solution(infile='ch4_plasma.cti', efile='lxcat.yaml')

    def test_electron_properties(self):
        self.gas.TPX = 1000, ct.one_atm, 'O2:1.0, E:1e-10'
        self.gas.electric_field = 1e5
        self.gas.electric_field_freq = 0.0
        self.assertNear(self.gas.electron_temperature, 15273, 1e-3)
        self.assertNear(self.gas.electron_mobility, 0.3981, 1e-4)
        self.assertNear(self.gas.electron_diffusivity, 0.5652, 1e-4)
        rate = self.gas.electron_rate_coefficient(19)
        self.assertNear(self.gas.net_plasma_production_rates[24] * 1e-10, rate, 1e-4)
        self.assertNear(self.gas.net_plasma_production_rates[24], 0.001877, 1e-3)
        self.assertNear(self.gas.electron_total_collision_frequency, 3.6957e11, 1e-3)
        self.assertNear(self.gas.electron_power_gain, 3.9811e9, 1e-3)
        self.assertNear(self.gas.electron_elastic_power_loss, 2.8671e7, 1e-3)
        self.assertNear(self.gas.mean_electron_energy, 1.9742, 1e-3)
        self.assertNear(self.gas.electric_field, 1e5, 1e-3)

    def test_change_electric_field_freq(self):
        self.gas.TPX = 1000, ct.one_atm, 'O2:1.0'
        self.gas.electric_field = 1e5
        self.gas.electric_field_freq = 0.0
        Te0 = self.gas.electron_temperature
        self.gas.electric_field_freq = 1e9
        Te = self.gas.electron_temperature
        self.gas.electric_field_freq = 2e9
        self.assertLess(Te, Te0)
        self.assertLess(self.gas.electron_temperature, Te)
        self.assertNear(Te, Te0, 1e-3)
        self.assertNear(self.gas.electron_temperature, Te, 1e-3)

    def test_change_gas_temperature(self):
        self.gas.TPX = 1000, ct.one_atm, 'O2:1.0'
        self.gas.electric_field = 1e5
        self.gas.electric_field_freq = 0.0
        Te0 = self.gas.electron_temperature
        self.gas.TP = 1100, ct.one_atm * 1.1
        # The gas temperature is important only when E/N is small
        self.assertNotEqual(self.gas.electron_temperature, Te0)
        self.assertNear(self.gas.electron_temperature, Te0, 1e-4)
        self.gas.electric_field = 1.0
        Te = self.gas.electron_temperature
        self.gas.TP = 1000, ct.one_atm
        self.assertLess(0.05, abs(self.gas.electron_temperature - Te) / Te)
        self.gas.electric_field = 1e5

    def test_change_electric_field_strength(self):
        self.gas.TPX = 1000, ct.one_atm, 'O2:1.0'
        self.gas.electric_field = 1e5
        self.gas.electric_field_freq = 0.0
        grid = np.linspace(0.0, 2.0, num=1000)
        self.gas.set_electron_energy_grid(grid)
        self.gas.electric_field = 1.0
        # The electron temperature approach gas temperature when E is small
        self.assertNear(self.gas.electron_temperature, self.gas.T, 1e-3)

