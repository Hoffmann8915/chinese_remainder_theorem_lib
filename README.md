# chinese_remainder_theorem_lib

Solve a system of simultaneous modular congruences with pairwise-coprime moduli using the extended Euclidean algorithm. Returns the unique solution modulo the product of the moduli.

## Usage

```python
from chinese_remainder_theorem_lib import solve

# Find x such that x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7)
x = solve([2, 3, 2], [3, 5, 7])
print(x)  # 23
```

`crt` is an alias for `solve`. `solve_all` returns the full solution plus each partial solution obtained by dropping one constraint.

## Exported names

| Name | Signature | Description |
|------|-----------|-------------|
| `solve(residues: list[int], moduli: list[int]) -> int` | `(residues, moduli)` | Main solver. Returns the unique solution in `[0, product)`. |
| `crt` | same as `solve` | Alias for `solve`. |
| `solve_all(residues: list[int], moduli: list[int]) -> list[int]` | `(residues, moduli)` | Full solution followed by each single-drop partial solution. |

## Why

Standard library `pow(a, -1, m)` exists in Python 3.8+ for modular inverse, but building the full CRT solution from scratch with an explicit extended Euclidean algorithm keeps the logic inspectable and dependency-free. The trade-off is that you must guarantee pairwise coprimality yourself — the library checks it and raises `ValueError` if violated, but it does not attempt to repair or decompose non-coprime systems.

## Edge cases

Negative residues are accepted and reduced modulo their modulus before solving. A modulus of 1 makes its congruence trivially satisfied (everything is 0 mod 1). Moduli must be positive integers; passing 0 or a negative value raises `ValueError`. The result is always in the range `[0, product_of_moduli)`.

## Performance

The window keeps a bounded buffer, so `push` is constant time and memory does not
grow with the length of the stream. `peak` and `trough` are linear in the window
size, which is the trade that keeps `push` cheap.

