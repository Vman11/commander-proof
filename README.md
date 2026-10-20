# GSR shuffles for a Commander deck — proof + code

**Question:** how well does riffle shuffling actually randomize a 100-card
Commander deck, and can that be *proved*, not just simulated?

**Answer:** yes — the Gilbert–Shannon–Reeds (GSR) model has an exact closed
form (Bayer & Diaconis, 1992), proved in full in
[`docs/PROOF.md`](docs/PROOF.md) and verified by brute force and simulation
in this repo.

## The result

For one $a$-shuffle of $n$ cards, the probability of any permutation $\pi$ is

$$Q_a(\pi) = \frac{1}{a^n}\binom{n + a - r(\pi)}{n},$$

depending on $\pi$ only through its number of rising sequences $r(\pi)$.
$k$ ordinary riffles are an $a = 2^k$-shuffle. Specializing to $n = 100$
and measuring total variation distance to uniform (computed **exactly**,
in rational arithmetic — `results/tv_commander.csv`):

| shuffles | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|
| TV distance to random | 0.747 | 0.429 | 0.224 | 0.112 | 0.057 | 0.028 |

**Bottom line:** the familiar "7 shuffles" rule is for a 52-card deck. A
Commander deck needs **~10 shuffles** to reach the cutoff
($\tfrac32\log_2 100 \approx 9.97$), and 11–12 to get the distance under
~0.05. The 99-card library (commander in the command zone) is essentially
identical.

## What's here

- `docs/PROOF.md` — the complete proof: physical definition → drop rule is
  a uniform interleaving (Lemma 1) → every (cut, interleaving) pair has
  probability exactly $a^{-n}$ (Lemma 2) → stars-and-bars count gives the
  formula (Theorem) → Commander corollaries, including the exact
  total-variation formula.
- `src/gsr.py` — exact probabilities as `Fraction`s, rising sequences, a
  literal physical simulator (binomial cut + size-proportional drops),
  Eulerian counts, and exact total variation for any $n$.
- `src/simulate_commander.py` — Monte-Carlo check at $n = 100$: simulated
  rising-sequence means vs. exact theory (results in
  `results/simulation.txt`; they agree to ~0.1).
- `tests/test_gsr.py` — brute-force verification: for $n \le 6$ and
  $a \in \{2, 3\}$, enumerate the physical process directly and compare
  with the formula on **every** permutation, exactly. Also reproduces the
  published $n = 52$ Bayer–Diaconis TV table as a control.
- `results/` — exact TV tables for $n = 100$ and $n = 99$.

## Run it

```bash
python3 tests/test_gsr.py          # all checks, ~1s, no dependencies
python3 src/gsr.py                 # exact TV table for n = 100
python3 src/simulate_commander.py  # simulation vs. theory at n = 100
```

Pure Python 3 standard library; nothing to install.

## References

Bayer & Diaconis, *Trailing the Dovetail Shuffle to its Lair*, Ann. Appl.
Prob. 2 (1992) 294–313; Mann (1995) for the elementary counting proof.
Full citations in `docs/PROOF.md`.
