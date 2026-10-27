"""Core implementation of the Chinese Remainder Theorem.

This module provides functions to solve a system of modular congruences
of the form x ≡ r_i (mod m_i) where the moduli m_i are pairwise coprime.
The solution is unique modulo the product of all moduli.

The implementation uses the extended Euclidean algorithm to compute
modular inverses. Only built-in integer arithmetic is used, which
means results are exact (no floating-point rounding).
"""

from __future__ import annotations

__all__ = ["solve", "crt", "solve_all"]


def extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    """Return (gcd, x, y) such that a*x + b*y == gcd(a, b).

    The Bézout identity is always satisfied for any pair of integers,
    including negative inputs. This is a requirement for the CRT solver
    when negative residues are given.
    """
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t
    return old_r, old_s, old_t


def mod_inverse(a: int, m: int) -> int:
    """Return the modular inverse of a modulo m.

    Raises ValueError when gcd(a, m) != 1, i.e. the inverse does not
    exist. We reduce a modulo m first so that the caller may pass any
    representative of the residue class (including negatives).
    """
    a_reduced = a % m
    if a_reduced < 0:  # defensive; Python's % already yields non-negative
        a_reduced += m
    g, x, _ = extended_gcd(a_reduced, m)
    if g != 1:
        raise ValueError(f"{a} has no inverse modulo {m} (gcd={g})")
    return x % m


def _check_pairwise_coprime(moduli: list[int]) -> None:
    """Raise ValueError if any pair of moduli shares a common factor.

    The CRT uniqueness guarantee requires pairwise coprimality. We verify
    this explicitly rather than relying on the caller, because a silent
    wrong answer is far worse than a loud error.
    """
    n = len(moduli)
    for i in range(n):
        for j in range(i + 1, n):
            g, _, _ = extended_gcd(moduli[i], moduli[j])
            if abs(g) != 1:
                raise ValueError(
                    f"moduli {moduli[i]} and {moduli[j]} are not coprime "
                    f"(gcd={g})"
                )


def solve(residues: list[int], moduli: list[int]) -> int:
    """Solve the system x ≡ r_i (mod m_i) for pairwise-coprime moduli.

    Parameters
    ----------
    residues : list of int
        The remainders r_i. Any integer is accepted (positive, negative,
        or zero) because the congruence x ≡ r (mod m) depends only on
        the residue class of r modulo m.
    moduli : list of int
        The moduli m_i. They must be pairwise coprime. Each modulus must
        be positive; a non-positive modulus is rejected because the
        modular ring Z/mZ is only standardly defined for m ≥ 1, and m=0
        would make every congruence trivially true or meaningless.

    Returns
    -------
    int
        The unique solution in the range [0, product_of_moduli). The
        lower bound is inclusive; the upper bound is exclusive, so the
        result is a canonical representative of the solution class.

    Raises
    ------
    ValueError
        If the lists have different lengths, are empty, contain a
        non-positive modulus, or the moduli are not pairwise coprime.
    """
    if len(residues) != len(moduli):
        raise ValueError(
            f"residues and moduli must have the same length, "
            f"got {len(residues)} and {len(moduli)}"
        )
    if len(moduli) == 0:
        raise ValueError("at least one congruence is required")
    for m in moduli:
        if m <= 0:
            raise ValueError(f"moduli must be positive, got {m}")
    _check_pairwise_coprime(moduli)

    total_product = 1
    for m in moduli:
        total_product *= m

    result = 0
    for r, m in zip(residues, moduli):
        partial = total_product // m
        inv = mod_inverse(partial, m)
        result += (r % m) * partial * inv

    return result % total_product


def solve_all(residues: list[int], moduli: list[int]) -> list[int]:
    """Return [solve(...)] and all solutions obtained by dropping each constraint in turn.

    This is occasionally useful when one residue is suspected to be wrong
    and you want to see the partial solutions. The full solution is
    always the first element of the returned list.
    """
    full = solve(residues, moduli)
    results = [full]
    for i in range(len(residues)):
        sub_res = residues[:i] + residues[i + 1:]
        sub_mod = moduli[:i] + moduli[i + 1:]
        results.append(solve(sub_res, sub_mod))
    return results


crt = solve
"""Alias for solve."""
