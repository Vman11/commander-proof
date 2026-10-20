# The Gilbert–Shannon–Reeds model for a Commander deck: a complete proof

**Claim proved.** For a deck of $n$ distinct cards, one $a$-shuffle in the
Gilbert–Shannon–Reeds (GSR) model produces a permutation $\pi$ with probability

$$
Q_a(\pi) \;=\; \frac{1}{a^{n}}\binom{n + a - r(\pi)}{\,n\,},
$$

where $r(\pi)$ is the number of *rising sequences* of $\pi$ and the binomial
coefficient is $0$ when $a < r(\pi)$. Everything below is proved from the
physical definition of the shuffle, for general $n$ and $a$; the Commander
case is the specialization $n = 100$ (Corollary 4), with the exact number of
shuffles a 100-card deck needs in Corollary 5.

This is the main theorem of Bayer & Diaconis (1992). The proof given here is
the elementary counting proof (in the form exposited by Mann, 1995), written
out in full, and every step is checked by brute force in
`tests/test_gsr.py` for small decks.

---

## 1. Setup

Label the cards $1, \dots, n$ in their initial order (card $1$ on top).
A shuffle outcome is the sequence $\pi = (\pi_1, \dots, \pi_n)$ of labels,
top to bottom, after shuffling. Write $\operatorname{pos}_\pi(k)$ for the
position of card $k$ in $\pi$.

**Definition (rising sequences).** Call $k \in \{1, \dots, n-1\}$ a *bad gap*
(for $\pi$) if $\operatorname{pos}_\pi(k+1) < \operatorname{pos}_\pi(k)$,
i.e. card $k+1$ sits above card $k$. Cutting the label sequence
$1, 2, \dots, n$ at every bad gap partitions it into maximal intervals of
consecutive labels that appear in $\pi$ in increasing order. These intervals
are the **rising sequences** of $\pi$, and $r(\pi)$ is their number:

$$
r(\pi) \;=\; 1 + \#\{\text{bad gaps}\} \;=\; 1 + \operatorname{des}(\pi^{-1}).
$$

*Example.* $\pi = 7\,1\,6\,2\,4\,8\,3\,5\,9$ has bad gaps at $k=3,5,6$,
so its rising sequences are $123,\ 45,\ 6,\ 789$ and $r(\pi) = 4$.

**Definition (GSR $a$-shuffle).**

1. **Cut.** Cut the deck into $a$ consecutive piles, of sizes
   $n_1, \dots, n_a \ge 0$ (empty piles allowed), with multinomial probability
   $\displaystyle \Pr(n_1,\dots,n_a) = \frac{n!}{n_1! \cdots n_a!}\cdot\frac{1}{a^n}$.
   For $a = 2$ this is the familiar $\binom{n}{k}/2^n$ binomial cut.
2. **Riffle.** Repeatedly drop the next card from one of the piles, choosing
   pile $i$ with probability proportional to its current size $A_i$,
   i.e. with probability $A_i / (A_1 + \cdots + A_a)$, until all piles are
   exhausted. Each pile's internal order is preserved.

An ordinary riffle shuffle is the case $a = 2$.

## 2. Two lemmas

**Lemma 1 (the drop rule is a uniform interleaving).**
*Fix the pile sizes $n_1, \dots, n_a$. Under the GSR drop rule, every
interleaving (merge) of the piles that preserves each pile's internal order
occurs with the same probability, $\prod_i n_i! \,/\, n!$.*

*Proof.* Fix one interleaving. Reading it from first drop to last, it drops
from pile $i$ exactly $n_i$ times, and at those moments the number of cards
remaining in pile $i$ runs through $n_i, n_i - 1, \dots, 1$ in that order
(regardless of when the drops happen). Meanwhile the total number of cards
remaining runs through $n, n-1, \dots, 1$. So the probability of this
interleaving is

$$
\prod_{\text{drops}} \frac{\text{size of chosen pile}}{\text{cards remaining}}
\;=\; \frac{\prod_i (n_i \cdot (n_i-1) \cdots 1)}{n \cdot (n-1) \cdots 1}
\;=\; \frac{\prod_i n_i!}{n!},
$$

the same for every interleaving. ∎

**Lemma 2 (every (cut, interleaving) pair has probability exactly $a^{-n}$).**

*Proof.* Multiply the probability of the cut (step 1 of the definition) by
the probability of the interleaving given that cut (Lemma 1):

$$
\frac{n!}{\prod_i n_i!}\cdot\frac{1}{a^n} \;\cdot\; \frac{\prod_i n_i!}{n!}
\;=\; \frac{1}{a^n}.
$$

The pile sizes cancel completely. ∎

Lemma 2 is the engine of the whole theory: to find $Q_a(\pi)$ we only have
to **count** the (cut, interleaving) pairs that produce $\pi$, and divide
by $a^n$.

## 3. Main theorem (Bayer–Diaconis, 1992)

**Theorem.** $\displaystyle Q_a(\pi) = \frac{1}{a^n}\binom{n + a - r(\pi)}{n}$,
with the right-hand side read as $0$ when $a < r(\pi)$.

*Proof.* By Lemma 2,

$$
Q_a(\pi) = \frac{N(\pi)}{a^n},
\qquad
N(\pi) = \#\{(\text{cut}, \text{interleaving}) \text{ pairs producing } \pi\}.
$$

**Step 1: pairs producing $\pi$ are pile assignments.**
A cut is the same data as a weakly increasing map
$p : \{1,\dots,n\} \to \{1,\dots,a\}$, where $p_k$ is the pile containing
card $k$: weak increase is exactly "the piles are consecutive intervals of
the deck", and values of $\{1,\dots,a\}$ skipped by $p$ are the empty piles.
Given a cut $p$ and a target $\pi$, there is **at most one** interleaving
that yields $\pi$ — the cards must drop in the order $\pi$ lists them — and
it exists iff $\pi$ respects the internal order of every pile. Hence pairs
producing $\pi$ are in bijection with the weakly increasing maps $p$ such
that $\pi$ respects every pile's internal order.

**Step 2: which pile assignments work.**
A pile is an interval of consecutive labels, so $\pi$ respects its order iff
$\pi$ respects the order of each *adjacent pair* $k, k+1$ inside it
(transitivity along the interval). An adjacent pair is respected iff $k$ is
not a bad gap. Therefore $p$ works iff

$$
\text{for every bad gap } k:\quad p_{k+1} > p_k
\qquad\text{(a strict increase of at least 1),}
$$

while at all other gaps $p_{k+1} \ge p_k$ is automatic. In words: **every
bad gap — every boundary between two rising sequences — must be separated
by at least one pile boundary**, and there are no other constraints.

**Step 3: counting (stars and bars).**
There are $r - 1$ bad gaps. Define $b_k$ = the number of bad gaps strictly
before label $k$ (so $0 \le b_1 \le \dots \le b_n = r-1$, incrementing
exactly at bad gaps), and set $s_k = p_k - b_k$. The strict increases at bad
gaps become weak increases after this shift, the weak increases elsewhere
stay weak, so $s$ is an arbitrary weakly increasing map
$\{1,\dots,n\} \to \{1, \dots, a - (r-1)\}$, and $p \mapsto s$ is a bijection.
If $a < r$ no such $p$ exists and $N(\pi) = 0$, as claimed. Otherwise, the
number of weakly increasing maps from an $n$-element chain to a chain of
$m = a - r + 1$ values is the stars-and-bars count

$$
N(\pi) = \binom{n + m - 1}{n} = \binom{n + a - r(\pi)}{n},
$$

independent of everything about $\pi$ except $r(\pi)$. Dividing by $a^n$
(Lemma 2) completes the proof. ∎

**Sanity checks** (all also verified exactly in code — in Python in
`tests/test_gsr.py`, and independently in Lean 4 in `lean/GsrVerify.lean`,
which checks the formula against brute force on *every* permutation for
$(n, a) \in \{(4,2), (6,2), (4,3), (5,3)\}$ and pins the exact $n = 100$
total-variation values of Corollary 5 into 0.001-wide rational intervals):

- Summing over $\pi$ grouped by $r$ gives $\sum_\pi Q_a(\pi) = 1$
  (checked for $n \le 7$, several $a$).
- $a = 1$: only $\pi = \mathrm{id}$ (with $r = 1$) gets probability
  $\binom{n}{n}/1 = 1$. No shuffle, no change. ✓
- $a \to \infty$: $\binom{n+a-r}{n}/a^n \to 1/n!$ for every $\pi$: a huge
  number of piles is a uniform shuffle. ✓

## 4. Corollaries for a Commander deck ($n = 100$)

**Corollary 1 (one shuffle barely shuffles).** With $a = 2$, the theorem
gives: $\pi$ is reachable in one GSR shuffle iff $r(\pi) \le 2$; the
identity has probability $(n+1)/2^n$ and every other reachable $\pi$ has
probability $1/2^n$. For $n = 100$: the deck comes back completely
unshuffled with probability $101/2^{100} \approx 7.97 \times 10^{-29}$ —
tiny, but the real point is the support: after one shuffle, **almost every
ordering of a Commander deck is impossible**, because a random permutation
of 100 cards has $r \approx 50.5$ on average, not $\le 2$.

**Corollary 2 ($k$ shuffles are a $2^k$-shuffle).** Bayer & Diaconis prove
$Q_a * Q_b = Q_{ab}$ (convolution on the symmetric group): an $a$-shuffle
followed by an independent $b$-shuffle is an $ab$-shuffle. (Proof idea:
run the second shuffle's piles *inside* the first shuffle's piles in the
inverse description — refining $a$ piles into $ab$ piles multiplies the
pile counts and, by Lemma 2, every refined pair still has equal
probability.) Hence $k$ ordinary riffle shuffles give

$$
Q^{*k}(\pi) \;=\; \frac{1}{2^{kn}}\binom{n + 2^{k} - r(\pi)}{n}.
$$

**Corollary 3 (possibility threshold).** A permutation with $r(\pi) = 100$
(the fully reversed deck) first becomes possible when $2^k \ge 100$, i.e.
after $k = \lceil \log_2 100 \rceil = 7$ shuffles. Possible $\ne$ likely —
see Corollary 5.

**Corollary 4 (exact distance to random).** Grouping the theorem by $r$ and
using that $\pi \mapsto \pi^{-1}$ is a bijection taking $r(\pi)$ to
$\operatorname{des}(\pi) + 1$, the number of permutations of $n$ with a
given $r$ is the Eulerian number $A(n, r)$. So the total variation distance
to the uniform distribution $U$ after $k$ shuffles is the *exact* rational
number

$$
\|Q^{*k} - U\|_{\mathrm{TV}}
\;=\; \tfrac12 \sum_{r=1}^{n} A(n, r)\,
\Bigl|\, \tbinom{n + 2^k - r}{n} / 2^{kn} \;-\; 1/n! \,\Bigr|,
$$

computed with exact big-integer arithmetic in `src/gsr.py` (no floating
point anywhere until the final display). As a control, the same code at
$n = 52$ reproduces the published Bayer–Diaconis table exactly
($k=5{:}\ 0.924$, $6{:}\ 0.614$, $7{:}\ 0.334$, $8{:}\ 0.167$, $9{:}\ 0.085$).

**Corollary 5 (how many shuffles a Commander deck actually needs).**
Exact values for $n = 100$ (`results/tv_commander.csv`):

| shuffles $k$ | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|
| TV distance | 0.9815 | 0.7473 | 0.4288 | 0.2239 | 0.1122 | 0.0566 | 0.0284 | 0.0142 |

($k \le 5$ rounds to $1.000000$; e.g. $k = 5$ is $1 - 10^{-10}$, roughly.)
Once past the threshold the distance halves with each extra shuffle,
pinpointing the famous cutoff at $\tfrac32 \log_2 n = \tfrac32 \log_2 100
\approx 9.97$ shuffles. Concretely:

- **7 shuffles — the "standard advice" for a 52-card deck — leave a
  Commander deck at TV distance 0.75: nowhere near random.**
- **10 shuffles** get to $\approx 0.11$, **11** to $\approx 0.057$,
  **12** to $\approx 0.028$.

Strictly, the shuffled object in a game is the 99-card library (the
commander starts in the command zone); the theorem is for general $n$ and
the $n = 99$ numbers in `results/tv_commander.csv` are nearly identical
(11 shuffles: $0.0556$).

## 5. Scope and honest limits

- This proves the distribution **of the GSR model**, i.e. of an idealized
  riffle whose cut is binomial and whose drops are size-proportional.
  Real thumbs are not exactly this; the model is the standard, empirically
  good description of riffle shuffling, not a theorem about your hands.
- "Random" here means uniform over all $100!$ orders, measured in total
  variation. Weaker notions (a single card's position, the opening seven)
  mix somewhat faster; TV is the conservative, standard choice.
- Mash and pile "shuffles" are different processes and are **not** covered
  by this theorem.

## References

- D. Bayer & P. Diaconis, *Trailing the Dovetail Shuffle to its Lair*,
  Annals of Applied Probability 2 (1992), 294–313.
  <https://www.stat.berkeley.edu/~aldous/157/Papers/bayer_diaconis.pdf>
- B. Mann, *How Many Times Should You Shuffle a Deck of Cards?* (1995) —
  elementary exposition of the counting proof above.
- E. Gilbert (1955), C. Shannon (unpublished memo), J. Reeds (1981) —
  independent origins of the model.
- Gilbert–Shannon–Reeds model, Wikipedia (overview and equivalent
  descriptions). <https://en.wikipedia.org/wiki/Gilbert%E2%80%93Shannon%E2%80%93Reeds_model>
