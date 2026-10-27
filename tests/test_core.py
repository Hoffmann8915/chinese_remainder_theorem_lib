import unittest

from chinese_remainder_theorem_lib import solve, crt, solve_all
from chinese_remainder_theorem_lib.core import extended_gcd, mod_inverse


class TestExtendedGCD(unittest.TestCase):
    """The extended Euclidean algorithm underpins the whole library,
    so we test it in isolation."""

    def test_gcd_positive(self):
        g, x, y = extended_gcd(30, 12)
        self.assertEqual(g, 6)
        self.assertEqual(30 * x + 12 * y, g)

    def test_gcd_coprime(self):
        g, x, y = extended_gcd(17, 13)
        self.assertEqual(g, 1)
        self.assertEqual(17 * x + 13 * y, g)

    def test_gcd_with_zero(self):
        g, x, y = extended_gcd(15, 0)
        self.assertEqual(g, 15)
        self.assertEqual(15 * x + 0 * y, g)

    def test_bezout_identity_holds(self):
        # The identity a*x + b*y = gcd must hold for all inputs.
        for a in range(-10, 11):
            for b in range(-10, 11):
                if a == 0 and b == 0:
                    continue
                g, x, y = extended_gcd(a, b)
                self.assertEqual(a * x + b * y, g)


class TestModInverse(unittest.TestCase):
    def test_basic_inverse(self):
        self.assertEqual((3 * mod_inverse(3, 7)) % 7, 1)

    def test_no_inverse(self):
        with self.assertRaises(ValueError):
            mod_inverse(2, 4)

    def test_negative_input(self):
        # -1 mod 7 is 6, whose inverse is 6 (6*6=36≡1 mod 7).
        inv = mod_inverse(-1, 7)
        self.assertEqual((-1 * inv) % 7, 1)


class TestSolve(unittest.TestCase):
    """Tests for the main CRT solver."""

    def test_classic_two_congruences(self):
        # x ≡ 2 (mod 3), x ≡ 3 (mod 5)  →  x = 23 (mod 15)
        result = solve([2, 3], [3, 5])
        self.assertEqual(result, 8)
        self.assertEqual(result % 3, 2)
        self.assertEqual(result % 5, 3)

    def test_three_congruences(self):
        # x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7) → 23 (mod 105)
        result = solve([2, 3, 2], [3, 5, 7])
        self.assertEqual(result, 23)
        self.assertEqual(result % 3, 2)
        self.assertEqual(result % 5, 3)
        self.assertEqual(result % 7, 2)

    def test_single_congruence(self):
        result = solve([4], [7])
        self.assertEqual(result, 4)

    def test_result_in_range(self):
        residues = [1, 2, 3, 4]
        moduli = [5, 7, 11, 13]
        result = solve(residues, moduli)
        product = 5 * 7 * 11 * 13
        self.assertGreaterEqual(result, 0)
        self.assertLess(result, product)

    def test_negative_residues(self):
        # x ≡ -1 (mod 5), x ≡ -1 (mod 7) → x = 34 (mod 35)
        result = solve([-1, -1], [5, 7])
        self.assertEqual(result, 34)
        self.assertEqual(result % 5, 4)
        self.assertEqual(result % 7, 6)

    def test_large_residue_wraps(self):
        # Residues larger than the modulus should be reduced.
        # 7 mod 5 = 2, 8 mod 7 = 1
        result = solve([7, 8], [5, 7])
        self.assertEqual(result % 5, 2)
        self.assertEqual(result % 7, 1)

    def test_zero_residues(self):
        result = solve([0, 0], [3, 5])
        self.assertEqual(result, 0)

    def test_result_is_zero_solution(self):
        # When all residues are 0 the answer is 0, not the product.
        result = solve([0, 0, 0], [3, 5, 7])
        self.assertEqual(result, 0)

    def test_alias_crt_exists(self):
        self.assertIs(crt, solve)

    def test_mismatched_lengths(self):
        with self.assertRaises(ValueError):
            solve([1, 2], [3, 5, 7])

    def test_empty_input(self):
        with self.assertRaises(ValueError):
            solve([], [])

    def test_non_positive_modulus(self):
        with self.assertRaises(ValueError):
            solve([1], [0])
        with self.assertRaises(ValueError):
            solve([1], [-3])

    def test_non_coprime_moduli(self):
        with self.assertRaises(ValueError):
            solve([1, 2], [4, 6])


class TestSolveAll(unittest.TestCase):
    def test_returns_full_plus_dropped(self):
        residues = [2, 3, 2]
        moduli = [3, 5, 7]
        results = solve_all(residues, moduli)
        self.assertEqual(len(results), 4)
        self.assertEqual(results[0], 23)
        # Dropping index 0: solve([3,2],[5,7])
        self.assertEqual(results[1] % 5, 3)
        self.assertEqual(results[1] % 7, 2)


class TestEdgeCases(unittest.TestCase):
    def test_modulus_one(self):
        # Every integer ≡ 0 (mod 1), so this constraint is trivial.
        # The answer is determined by the remaining congruence.
        result = solve([0, 3], [1, 5])
        self.assertEqual(result, 3)

    def test_large_numbers(self):
        residues = [1, 2, 3]
        moduli = [9999991, 99999931, 99999989]
        result = solve(residues, moduli)
        self.assertEqual(result % 9999991, 1)
        self.assertEqual(result % 99999931, 2)
        self.assertEqual(result % 99999989, 3)


if __name__ == "__main__":
    unittest.main()
