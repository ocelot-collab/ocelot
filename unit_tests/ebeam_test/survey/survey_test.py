import unittest
import numpy as np
import sys
import os

# --- IMPORT SETUP ---
# Adjust this path if necessary so Python can find the 'ocelot' package.
# If you run this file from the project root, this usually isn't needed.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import your Ocelot classes
# (Adjust these import paths to match exactly where your files are located)
from ocelot.cpbd.elements import Drift, SBend, Hcor
from ocelot.cpbd.magnetic_lattice import MagneticLattice


class TestSurveyDrift(unittest.TestCase):
    def test_straight_line(self):
        """Checks that Drifts propagate purely along Z/S with no deviation."""
        drifts = [Drift(l=1.0, eid=f"d{i}") for i in range(10)]
        lat = MagneticLattice(drifts)

        mid_data, end_data = lat.survey(X0=0, Y0=0, Z0=0)
        last_point = end_data[-1]

        self.assertAlmostEqual(last_point['X'], 0.0)
        self.assertAlmostEqual(last_point['Y'], 0.0)
        self.assertAlmostEqual(last_point['Z'], 10.0)
        self.assertAlmostEqual(last_point['S'], 10.0)
        self.assertAlmostEqual(last_point['THETA'], 0.0)
        self.assertAlmostEqual(last_point['PHI'], 0.0)


class TestSurveyBend(unittest.TestCase):
    def test_90_deg_bend(self):
        """Checks a standard horizontal 90-degree bend."""
        L = np.pi / 2.0
        angle = np.pi / 2.0
        bend = SBend(l=L, angle=angle)
        lat = MagneticLattice([bend])

        mid_data, end_data = lat.survey()

        # 1. End Point Check (Last item in end_data)
        last = end_data[-1]
        self.assertAlmostEqual(last['X'], -1.0)  # Bend Right -> Negative X
        self.assertAlmostEqual(last['Y'], 0.0)
        self.assertAlmostEqual(last['Z'], 1.0)

        # Angle should be -pi/2 (Right turn)
        self.assertAlmostEqual(last['THETA'], -np.pi / 2.0)

        # 2. Midpoint Check
        # FIX: Access [-1] (The Bend), NOT [0] (The Start Point)
        mid = mid_data[-1]

        expected_x = 1.0 * (np.cos(np.pi / 4) - 1.0)  # approx -0.29289
        expected_z = 1.0 * np.sin(np.pi / 4)  # approx 0.7071

        self.assertAlmostEqual(mid['X'], expected_x)
        self.assertAlmostEqual(mid['Z'], expected_z)


class TestSurveyVertical(unittest.TestCase):
    def test_vertical_bend(self):
        """Checks a bend rotated 90 degrees via TILT."""
        L = np.pi / 2.0
        angle = np.pi / 2.0
        # TILT = 90 degrees -> Vertical Bend
        v_bend = SBend(l=L, angle=angle, tilt=np.pi / 2.0)

        lat = MagneticLattice([v_bend])
        _, end_data = lat.survey()
        last = end_data[-1]

        # X should stay 0, Y should take the displacement
        self.assertAlmostEqual(last['X'], 0.0)
        self.assertAlmostEqual(last['Y'], -1.0)
        self.assertAlmostEqual(last['Z'], 1.0)

        # CORRECTED Check:
        # 1. Phi (Elevation) must be -90 degrees (Straight Down)
        self.assertAlmostEqual(last['PHI'], -np.pi / 2.0)

        # 2. Theta (Azimuth) is undefined at the pole (0/0).
        # Instead of checking Theta, check the Direction Vectors (XPD/ZPD)
        # XPD and ZPD should both be 0 for a purely vertical bend.
        self.assertAlmostEqual(last['XPD'], 0.0)
        self.assertAlmostEqual(last['ZPD'], 0.0)


class TestSurveyCorrector(unittest.TestCase):
    def test_corrector_geometry(self):
        """Checks that Correctors have straight geometry despite having an angle."""
        cor = Hcor(l=1.0, angle=0.5)
        lat = MagneticLattice([cor])
        _, end_data = lat.survey()
        last = end_data[-1]

        # Geometry must be straight (Drift-like)
        self.assertAlmostEqual(last['X'], 0.0)
        self.assertAlmostEqual(last['Y'], 0.0)
        self.assertAlmostEqual(last['Z'], 1.0)
        self.assertAlmostEqual(last['THETA'], 0.0)


class TestSurveyLongList(unittest.TestCase):
    def test_adapter_logic(self):
        """Checks the translation from MAD8 physics to Engineering LongList format."""
        # Horizontal bend
        bend = SBend(l=1.0, angle=0.1)
        lat = MagneticLattice([bend])

        mid_std, _ = lat.survey()
        mid_long, _ = lat.survey_longlist(X0=0, Y0=0, Z0=0, theta0=0, phi0=0, chi0=0)

        item_std = mid_std[0]
        item_long = mid_long[0]

        # PSI renamed to CHI
        self.assertNotIn('PSI', item_long)
        self.assertIn('CHI', item_long)
        self.assertEqual(item_long['CHI'], item_std['PSI'])

        # Angles Swapped/Inverted
        self.assertEqual(item_long['THETA'], item_std['PHI'])
        self.assertEqual(item_long['PHI'], -1 * item_std['THETA'])

        # Cosines Swapped
        self.assertEqual(item_long['XPD'], item_std['YPD'])
        self.assertEqual(item_long['YPD'], -1 * item_std['XPD'])


class TestSurveyRing(unittest.TestCase):
    def test_closed_ring(self):
        """Checks that a full circle returns to 0,0,0."""
        # 4 Bends, 90 deg each
        bend = SBend(l=np.pi / 2, angle=np.pi / 2)
        ring_seq = [bend, bend, bend, bend]

        lat = MagneticLattice(ring_seq)
        _, end_data = lat.survey()
        last = end_data[-1]

        self.assertAlmostEqual(last['X'], 0.0)
        self.assertAlmostEqual(last['Y'], 0.0)
        self.assertAlmostEqual(last['Z'], 0.0)

        # Check direction is back to Z-axis
        self.assertAlmostEqual(last['XPD'], 0.0)
        self.assertAlmostEqual(last['ZPD'], 1.0)


if __name__ == '__main__':
    # This line runs all methods starting with 'test_' in all classes above
    unittest.main()
