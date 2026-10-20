"""
Gilbert-Shannon-Reeds (GSR) riffle shuffles: exact probabilities and simulation.

Core result implemented here (Bayer & Diaconis, 1992):
    For an a-shuffle of n distinct cards, the probability of a permutation pi is

        Q_a(pi) = C(n + a - r(pi), n) / a**n

    where r(pi) is the number of rising sequences of pi, and the binomial
    coefficient is taken to be 0 when a < r(pi).

    k ordinary GSR (2-) shuffles are a single a-shuffle with a = 2**k.

A Commander deck in Magic: The Gathering has n = 100 cards
(99-card library + 1 commander).  All exact routines below use Python
integers / Fractions, so the n = 100 values are exact rational numbers,
not floating-point approximations.
"""

from __future__ import annotations

import random
from fractions import Fraction
from math import comb, factorial


# ---------------------------------------------------------------------------
# Rising sequences
# ---------------------------------------------------------------------------

def rising_sequences(perm) -> int:
    """Number of rising sequences r(pi) of a permutation of 1..n.

    r(pi) = 1 + #{k in 1..n-1 : position(k+1) < position(k) in pi},
    i.e. 1 + the number of descents of pi^{-1}.
    """
    pos = {card: i for i, card in enumerate(perm)}
    n = len(perm)
    return 1 + sum(1 for k in range(1, n) if pos[k + 1] < pos[k])


def rising_sequence_blocks(perm) -> list[list[int]]:
    """The rising sequences themselves: maximal intervals of consecutive
    labels that appear in increasing order in `perm`."""
    pos = {card: i for i, card in enumerate(perm)}
    blocks, current = [], [1]
    for k in range(1, len(perm)):
        if pos[k + 1] < pos[k]:
            blocks.append(current)
            current = [k + 1]
        else:
            current.append(k + 1)
    blocks.append(current)
    return blocks


# ---------------------------------------------------------------------------
# Exact GSR probabilities
# ---------------------------------------------------------------------------

def q_a(perm, a: int) -> Fraction:
    """Exact probability of `perm` after one a-shuffle (Bayer-Diaconis)."""
    n = len(perm)
    r = rising_sequences(perm)
    if a < r:
        return Fraction(0)
    return Fraction(comb(n + a - r, n), a ** n)


def q_after_k_shuffles(perm, k: int) -> Fraction:
    """Exact probability of `perm` after k ordinary GSR 2-shuffles."""
    return q_a(perm, 2 ** k)


# ---------------------------------------------------------------------------
# Simulation of the physical GSR process
# ---------------------------------------------------------------------------

def a_shuffle(deck, a: int, rng: random.Random) -> list:
    """One a-shuffle of `deck`, simulated by the physical definition:

    1. Cut into a piles with multinomial(1/a, ..., 1/a) pile sizes.
    2. Drop cards from the piles one at a time, choosing a pile with
       probability proportional to its current size (the GSR drop rule).

    By Lemma 1 in docs/PROOF.md, step 2 is the same as picking a uniformly
    random interleaving of the piles; this function follows the drop rule
    literally so the simulation matches the physical description.
    """
    n = len(deck)
    # Multinomial pile sizes: assign each cut point independently.
    sizes = [0] * a
    for _ in range(n):
        sizes[rng.randrange(a)] += 1
    # Note: assigning pile *sizes* this way is NOT assigning cards to piles;
    # the piles are consecutive blocks of the deck, of the sizes drawn above.
    piles, i = [], 0
    for s in sizes:
        piles.append(list(deck[i:i + s]))
        i += s
    out = []
    piles = [p for p in piles]  # keep empties; the drop rule skips them
    remaining = sum(len(p) for p in piles)
    while remaining:
        # Choose a pile with probability proportional to its current size.
        x = rng.randrange(remaining)
        acc = 0
        for p in piles:
            if len(p) > x - acc:
                out.append(p.pop(0))
                break
            acc += len(p)
        remaining -= 1
    return out


def riffle_shuffle(deck, rng: random.Random) -> list:
    """One ordinary GSR shuffle (a = 2)."""
    return a_shuffle(deck, 2, rng)


def shuffle_n(deck, k: int, rng: random.Random) -> list:
    """k ordinary GSR shuffles."""
    for _ in range(k):
        deck = riffle_shuffle(deck, rng)
    return deck


# ---------------------------------------------------------------------------
# Exact total variation distance to uniform, for any n (incl. n = 100)
# ---------------------------------------------------------------------------

def eulerian_distribution(n: int) -> list[int]:
    """A(n, r) = number of permutations of n with exactly r rising sequences,
    for r = 1..n (index 0 of the returned list is unused / 0).

    Rising sequences of pi are descents of pi^{-1} (+1), and inversion is a
    bijection on permutations, so A(n, r) is the classical Eulerian number
    of permutations of n with r-1 descents.  Computed by the recurrence
        A(n, r) = r * A(n-1, r) + (n - r + 1) * A(n-1, r-1).
    Exact big-integer arithmetic throughout.
    """
    prev = [0] * (n + 2)
    prev[1] = 1  # n = 1
    for m in range(2, n + 1):
        cur = [0] * (n + 2)
        for r in range(1, m + 1):
            cur[r] = r * prev[r] + (m - r + 1) * prev[r - 1]
        prev = cur
    return prev


def total_variation(n: int, a: int) -> Fraction:
    """Exact total variation distance between Q_a and uniform on S_n.

    Groups permutations by r = r(pi): there are A(n, r) of them, each with
    probability C(n + a - r, n) / a**n (or 0 if a < r), vs. 1/n! uniform.
        TV = (1/2) * sum_r A(n, r) * | C(n + a - r, n)/a**n - 1/n! |
    """
    counts = eulerian_distribution(n)
    uniform = Fraction(1, factorial(n))
    total = Fraction(0)
    for r in range(1, n + 1):
        if counts[r] == 0:
            continue
        p = Fraction(comb(n + a - r, n), a ** n) if a >= r else Fraction(0)
        total += counts[r] * abs(p - uniform)
    return total / 2


def total_variation_after_k_shuffles(n: int, k: int) -> Fraction:
    """Exact TV distance to uniform after k ordinary GSR shuffles."""
    return total_variation(n, 2 ** k)


if __name__ == "__main__":
    n = 100  # a Commander deck
    print(f"Exact total variation distance to uniform, n = {n} (Commander)")
    for k in range(0, 16):
        tv = total_variation_after_k_shuffles(n, k)
        print(f"  k = {k:2d} shuffles:  TV = {float(tv):.6f}")
