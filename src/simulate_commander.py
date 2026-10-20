"""Simulation check for n = 100: physical GSR shuffles vs. exact theory.

For each number k of riffle shuffles, we simulate the physical process
(cut + size-proportional drops, repeated k times) and compare the empirical
mean number of rising sequences with the exact expectation implied by the
Bayer-Diaconis formula:

    E[r] = sum_r  r * A(100, r) * C(100 + 2**k - r, 100) / 2**(100k)

Run:  python3 src/simulate_commander.py [trials]
"""

import random
import sys
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gsr import eulerian_distribution, rising_sequences, shuffle_n  # noqa: E402

N = 100


def exact_expected_rising_sequences(n: int, k: int) -> float:
    a = 2 ** k
    counts = eulerian_distribution(n)
    total = 0.0
    for r in range(1, n + 1):
        if counts[r] and a >= r:
            total += r * counts[r] * comb(n + a - r, n) / a ** n
    return total


def main(trials: int = 300) -> None:
    rng = random.Random(20261020)
    print(f"n = {N} (Commander deck), {trials} simulated decks per k")
    print(f"{'k':>3} {'exact E[r]':>12} {'simulated E[r]':>15}")
    for k in (1, 3, 5, 7, 9, 11):
        exact = exact_expected_rising_sequences(N, k)
        sim = sum(
            rising_sequences(shuffle_n(list(range(1, N + 1)), k, rng))
            for _ in range(trials)
        ) / trials
        print(f"{k:>3} {exact:>12.2f} {sim:>15.2f}")
    print("\nUniform-random benchmark: E[r] = (n + 1) / 2 = 50.50")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 300)
