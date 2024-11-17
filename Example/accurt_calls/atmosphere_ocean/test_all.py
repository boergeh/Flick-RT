import unittest

class test_all(unittest.TestCase):
    def test_run(self):
        try:
            import compute_radiation
            import plot_radiation
            import plot_radiation_ratio
        except Exception as e:
            self.fail(f"Flick exception: {e}")

if __name__ == "__main__":
    unittest.main()
