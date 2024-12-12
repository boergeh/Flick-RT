import unittest

class test_all_python_plots(unittest.TestCase):
    def test_run(self):
        try:
            import cie_chromaticity
            import colorimitry
            import marine_iops
            import water_and_ice_absorption
            import material_angular_scattering
            import mie_asymmetry_factor
            import mie_scattering_matrix
            import optical_thickness
            import photon_counting
            import snow_asymmetry_factor
            import snow_ssalb
        except Exception as e:
            self.fail(f"Example python plot exception: {e}")

if __name__ == "__main__":
    unittest.main()
