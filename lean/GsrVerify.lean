/-!
# GSR shuffles for a Commander deck — machine-checked verification in Lean 4

This file verifies, by computation inside Lean, the quantitative claims of
`docs/PROOF.md` (the Bayer–Diaconis formula for the Gilbert–Shannon–Reeds
model, specialized to a 100-card Commander deck):

1. **Formula vs. pile assignments** (`formulaCheck`): for every permutation
   of small decks, the number of (cut) pile assignments compatible with the
   permutation — counted by brute force — equals the formula's count
   `choose (n + a - r) n`, where `r` is the number of rising sequences.
   Since every (cut, interleaving) pair has probability exactly `a⁻ⁿ`
   (Lemma 2 of the proof), equality of counts *is* equality of the
   distributions.  Checked for all permutations at
   (n, a) = (4, 2), (6, 2), (4, 3), (5, 3).
2. **Formula vs. the physical process** (`physicalCheck`): independent
   enumeration of cuts and interleavings (the literal cut-and-riffle
   process) gives the same count, for n = 4 (kernel `decide`) and n = 5
   (`native_decide`).
3. **The counts form a distribution** (`sumCheck`): the formula's counts
   over all permutations sum to exactly `aⁿ`.
4. **Commander headline numbers** (`tv100_*`): the exact total variation
   distance to uniform at n = 100, computed in rational arithmetic from the
   formula and Eulerian counts, lies in the 0.001-wide intervals matching
   `results/tv_commander.csv`, for k = 7 … 12 shuffles (a = 2ᵏ).

Scope, honestly: these are *computational* verifications (the large ones
via `native_decide`, i.e. verified modulo Lean's native-code reduction
axiom `Lean.ofReduceBool`; the small ones are pure kernel `decide`).
The general-n analytic proof is the paper proof in `docs/PROOF.md`;
formalizing that counting argument symbolically for all n in Lean is
future work and is *not* claimed here.

No Mathlib: only Lean core is used, so `lean GsrVerify.lean` suffices.
-/

/-! ## Small utilities -/

/-- List lookup with default 0 (indices here are always in range). -/
def nthOr : List Nat → Nat → Nat
  | [], _ => 0
  | x :: _, 0 => x
  | _ :: xs, i + 1 => nthOr xs i

/-- Factorial (Lean core no longer ships `Nat.factorial`, so we define it). -/
def fact : Nat → Nat
  | 0 => 1
  | n + 1 => (n + 1) * fact n

/-- Binomial coefficient `C(n, k)`, computed by the multiplicative formula
    with exact division at every step (each partial value is itself a
    binomial coefficient, hence an integer). `0` when `k > n`. -/
def choose (n k : Nat) : Nat :=
  if k > n then 0
  else
    let k' := min k (n - k)
    (List.range k').foldl (fun acc i => acc * (n - k' + 1 + i) / (i + 1)) 1

/-- The deck `[1, 2, …, n]`. -/
def deckOf (n : Nat) : List Nat := (List.range n).map (· + 1)

/-- 0-based position of card `x` in `perm`. -/
def posOf (x : Nat) : List Nat → Nat
  | [] => 0
  | y :: ys => if x = y then 0 else posOf x ys + 1

/-! ## Rising sequences -/

/-- The number of rising sequences `r(π) = 1 + #{k < n : pos(k+1) < pos(k)}`. -/
def rising (perm : List Nat) : Nat :=
  1 + ((List.range (perm.length - 1)).filter
        (fun i => posOf (i + 2) perm < posOf (i + 1) perm)).length

/-! ## All permutations -/

def insertAll (x : Nat) : List Nat → List (List Nat)
  | [] => [[x]]
  | y :: ys => (x :: y :: ys) :: (insertAll x ys).map (fun zs => y :: zs)

def permsOf : List Nat → List (List Nat)
  | [] => [[]]
  | x :: xs => (permsOf xs).flatMap (insertAll x)

/-! ## Pile assignments (cuts) and the count in the formula -/

/-- All weakly increasing lists of length `n` with values in `[lo, a)`:
    exactly the cuts of an n-card deck into `a` consecutive (possibly
    empty) piles, where entry `i` is the pile of card `i + 1`. -/
def weakIncFrom : Nat → Nat → Nat → List (List Nat)
  | 0, _, _ => [[]]
  | n + 1, lo, a =>
      ((List.range a).filter (fun v => lo ≤ v)).flatMap
        (fun v => (weakIncFrom n v a).map (fun rest => v :: rest))

/-- Is a list of distinct labels strictly increasing? -/
def increasing : List Nat → Bool
  | [] => true
  | [_] => true
  | x :: y :: ys => (x < y) && increasing (y :: ys)

/-- A pile assignment is *compatible* with `perm` if, for every pile, the
    cards of that pile appear in `perm` in increasing order — precisely
    the condition that some interleaving of the piles produces `perm`. -/
def compat (p perm : List Nat) (a : Nat) : Bool :=
  (List.range a).all
    (fun v => increasing (perm.filter (fun card => nthOr p (card - 1) == v)))

/-- Brute-force count: the number of cuts compatible with `perm`. -/
def bruteCount (n a : Nat) (perm : List Nat) : Nat :=
  ((weakIncFrom n 0 a).filter (fun p => compat p perm a)).length

/-- The Bayer–Diaconis count: `C(n + a - r(π), n)`. -/
def formulaCount (n a : Nat) (perm : List Nat) : Nat :=
  choose (n + a - rising perm) n

/-! ## The literal physical process: cuts + interleavings -/

/-- All interleavings (order-preserving merges) of two piles.
    Written with explicit fuel (structural recursion on a natural number)
    rather than well-founded recursion, so that the kernel itself — not
    just native code — can evaluate it in the `decide` checks below. -/
def mergesFuel : Nat → List Nat → List Nat → List (List Nat)
  | 0, _, _ => []
  | _, [], ys => [ys]
  | _, xs, [] => [xs]
  | f + 1, x :: xs, y :: ys =>
      (mergesFuel f xs (y :: ys)).map (x :: ·) ++ (mergesFuel f (x :: xs) ys).map (y :: ·)

def merges (xs ys : List Nat) : List (List Nat) :=
  mergesFuel (xs.length + ys.length) xs ys

/-- Number of (cut, interleaving) pairs of the physical 2-shuffle process
    that produce `perm`: sum over all binomial cuts `k = 0 … n`. -/
def physicalCount (n : Nat) (perm : List Nat) : Nat :=
  ((List.range (n + 1)).map (fun k =>
    (merges (deckOf k) ((List.range (n - k)).map (fun i => k + 1 + i))).countP
      (fun q => q == perm))).foldl (· + ·) 0

/-! ## Checks: formula vs. brute force -/

def formulaCheck (n a : Nat) : Bool :=
  (permsOf (deckOf n)).all (fun perm => bruteCount n a perm == formulaCount n a perm)

def physicalCheck (n : Nat) : Bool :=
  (permsOf (deckOf n)).all (fun perm => physicalCount n perm == formulaCount n 2 perm)

def sumCheck (n a : Nat) : Bool :=
  ((permsOf (deckOf n)).map (formulaCount n a)).foldl (· + ·) 0 == a ^ n

-- Pure kernel checks (no native code involved).
theorem formula_n4_a2_kernel : formulaCheck 4 2 = true := by decide
theorem formula_n4_a3_kernel : formulaCheck 4 3 = true := by decide
theorem physical_n4_kernel : physicalCheck 4 = true := by decide
theorem sum_n4_a2_kernel : sumCheck 4 2 = true := by decide

-- Larger exhaustive checks, via native-code evaluation.
theorem formula_n6_a2 : formulaCheck 6 2 = true := by native_decide
theorem formula_n5_a3 : formulaCheck 5 3 = true := by native_decide
theorem physical_n5 : physicalCheck 5 = true := by native_decide
theorem sum_n6_a2 : sumCheck 6 2 = true := by native_decide
theorem sum_n5_a3 : sumCheck 5 3 = true := by native_decide

/-! ## Exact total variation distance at n = 100 -/

/-- One row of Eulerian numbers: from row `m - 1` (= `prev`, entry `r` is
    `A(m-1, r)`, entry 0 unused) build row `m` by
    `A(m, r) = r·A(m-1, r) + (m - r + 1)·A(m-1, r-1)`. -/
def eulerianStep (m : Nat) (prev : List Nat) : List Nat :=
  (List.range (m + 1)).map
    (fun r => if r = 0 then 0
              else r * nthOr prev r + (m + 1 - r) * nthOr prev (r - 1))

def eulerianIter : Nat → List Nat → List Nat
  | 0, row => row
  | s + 1, row => eulerianIter s (eulerianStep row.length row)

/-- `eulerianCounts n`, entry `r`, is the Eulerian number `A(n, r)`: the
    number of permutations of `n` with exactly `r` rising sequences. -/
def eulerianCounts (n : Nat) : List Nat :=
  if n = 0 then [1] else eulerianIter (n - 1) [0, 1]

def tvTerm (counts : List Nat) (n a r : Nat) (uni : Rat) : Rat :=
  let c : Rat := (nthOr counts r : Rat)
  let p : Rat := if r ≤ a
    then (choose (n + a - r) n : Rat) / ((a ^ n : Nat) : Rat)
    else 0
  c * (if p ≥ uni then p - uni else uni - p)

/-- Exact total variation distance between the `a`-shuffle distribution
    and the uniform distribution on permutations of `n` cards:
    `½ · Σᵣ A(n, r) · |C(n + a - r, n)/aⁿ - 1/n!|`. -/
def tvRat (n a : Nat) : Rat :=
  let counts := eulerianCounts n
  let uni : Rat := 1 / (fact n : Rat)
  ((List.range (n + 1)).foldl (fun acc r => acc + tvTerm counts n a r uni) 0) / 2

-- Control at n = 52 (published Bayer–Diaconis values, k = 7 → ≈ 0.334).
theorem tv52_k7 : (334 / 1000 : Rat) ≤ tvRat 52 128 ∧ tvRat 52 128 < (335 / 1000 : Rat) := by
  native_decide

-- The Commander numbers (n = 100), matching results/tv_commander.csv.
theorem tv100_k7 : (747 / 1000 : Rat) ≤ tvRat 100 128 ∧ tvRat 100 128 < (748 / 1000 : Rat) := by
  native_decide
theorem tv100_k8 : (428 / 1000 : Rat) ≤ tvRat 100 256 ∧ tvRat 100 256 < (429 / 1000 : Rat) := by
  native_decide
theorem tv100_k9 : (223 / 1000 : Rat) ≤ tvRat 100 512 ∧ tvRat 100 512 < (224 / 1000 : Rat) := by
  native_decide
theorem tv100_k10 : (112 / 1000 : Rat) ≤ tvRat 100 1024 ∧ tvRat 100 1024 < (113 / 1000 : Rat) := by
  native_decide
theorem tv100_k11 : (56 / 1000 : Rat) ≤ tvRat 100 2048 ∧ tvRat 100 2048 < (57 / 1000 : Rat) := by
  native_decide
theorem tv100_k12 : (28 / 1000 : Rat) ≤ tvRat 100 4096 ∧ tvRat 100 4096 < (29 / 1000 : Rat) := by
  native_decide

-- What the computational theorems rest on (expected: `Lean.ofReduceBool`
-- for the `native_decide` ones, no axioms beyond the usual for `decide`).
#print axioms tv100_k10
#print axioms formula_n4_a2_kernel
