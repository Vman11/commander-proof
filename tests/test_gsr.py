"""Verification of the Bayer-Diaconis formula against brute force.

The brute-force side knows nothing about rising sequences: it enumerates
the physical GSR process directly (every cut into consecutive piles, and
every interleaving of the piles), weighting each (cut, interleaving) pair
by its true probability.  The formula side is `gsr.q_a`.  They must agree
exactly, as Fractions, on every permutation for the tested (n, a).

Run with:  python3 -m pytest tests/test_gsr.py   (or plain python3)
"""

import itertools
import random
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from gsr import (  # noqa: E402
    a_shuffle,
    eulerian_distribution,
    q_a,
    rising_sequences,
    total_variation,
)


def brute_force_physical_distribution(n: int, a: int) -> dict:
    """Exact distribution of one a-shuffle, by enumerating the process.

    A cut is a weakly increasing assignment p_1 <= ... <= p_n of the cards
    (in their original order) to the a consecutive piles; empty piles are
    skipped pile numbers.  Its sizes occur with multinomial probability
    n! / prod(n_i!) / a**n.  Given the sizes, the GSR drop rule picks each
    interleaving with probability prod(n_i!) / n! (Lemma 1 of the proof).
    Hence every (cut, interleaving) pair has probability exactly 1 / a**n.
    """
    dist = defaultdict(Fraction)
    for piles in itertools.combinations_with_replacement(range(a), n):
        for perm in itertools.permutations(range(1, n + 1)):
            pos = {c: i for i, c in enumerate(perm)}
            last, ok = {}, True
            for card, pile in zip(range(1, n + 1), piles):
                if pile in last and pos[card] < last[pile]:
                    ok = False
                    break
                last[pile] = pos[card]
            if ok:
                dist[perm] += Fraction(1, a ** n)
    return dist


def test_formula_matches_brute_force():
    for n in range(2, 7):
        for a in (2, 3):
            brute = brute_force_physical_distribution(n, a)
            for perm in itertools.permutations(range(1, n + 1)):
                assert q_a(perm, a) == brute.get(perm, Fraction(0)), (n, a, perm)
            assert sum(brute.values()) == 1


def test_formula_is_a_probability_distribution():
    for n in range(2, 8):
        for a in (2, 3, 4, 8):
            total = sum(
                q_a(p, a) for p in itertools.permutations(range(1, n + 1))
            )
            assert total == 1


def test_one_shuffle_has_at_most_two_rising_sequences():
    rng = random.Random(20261020)
    for _ in range(2000):
        out = a_shuffle(list(range(1, 101)), 2, rng)
        assert sorted(out) == list(range(1, 101))
        assert rising_sequences(out) <= 2


def test_eulerian_distribution_sums_to_factorial():
    from math import factorial

    for n in range(1, 13):
        counts = eulerian_distribution(n)
        assert sum(counts[1:]) == factorial(n)


def test_total_variation_52_matches_bayer_diaconis_table():
    # Published values (Bayer & Diaconis 1992), to 3 decimals.
    known = {5: 0.924, 6: 0.614, 7: 0.334, 8: 0.167, 9: 0.085}
    for k, tv in known.items():
        assert abs(float(total_variation(52, 2 ** k)) - tv) < 0.001


def test_simulation_matches_formula_small_deck():
    # Statistical check: n = 6, one shuffle, 60k trials, seeded.
    # Permutations split into r = 1 (prob 7/64), r = 2 (prob 1/64 each),
    # r >= 3 (prob 0).  Compare the empirical frequency of the identity
    # (the unique r = 1 permutation) with 7/64 = 0.109375.
    rng = random.Random(7)
    trials = 60_000
    hits = sum(
        1 for _ in range(trials)
        if a_shuffle(list(range(1, 7)), 2, rng) == list(range(1, 7))
    )
    assert abs(hits / trials - 7 / 64) < 0.01


if __name__ == "__main__":
    test_formula_matches_brute_force()
    test_formula_is_a_probability_distribution()
    test_one_shuffle_has_at_most_two_rising_sequences()
    test_eulerian_distribution_sums_to_factorial()
    test_total_variation_52_matches_bayer_diaconis_table()
    test_simulation_matches_formula_small_deck()
    print("All tests passed.")
