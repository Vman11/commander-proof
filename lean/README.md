# Lean 4 verification

`GsrVerify.lean` machine-checks the quantitative claims of this repo by
computation inside Lean 4 (core only — no Mathlib, no dependencies).

## Run it

```bash
# any Lean 4 toolchain (developed against 4.35.0)
lean GsrVerify.lean        # ~5 s; success is silence except the two
                           # `#print axioms` reports at the end
```

## What is checked

| Theorem(s) | Claim | Method |
|---|---|---|
| `formula_n4_a2_kernel`, `formula_n4_a3_kernel` | For **every** permutation of 4 cards (a = 2 and a = 3), the brute-force count of compatible pile assignments equals `C(n + a − r(π), n)` | kernel `decide` — depends on **no axioms** |
| `formula_n6_a2`, `formula_n5_a3` | Same statement, exhaustively, at n = 6 / n = 5 | `native_decide` |
| `physical_n4_kernel`, `physical_n5` | Same counts obtained by independently enumerating the literal physical process (every binomial cut × every interleaving) | `decide` / `native_decide` |
| `sum_n4_a2_kernel`, `sum_n6_a2`, `sum_n5_a3` | The formula's counts over all permutations sum to exactly `aⁿ` (it is a probability distribution) | `decide` / `native_decide` |
| `tv52_k7` | Control: exact TV distance at n = 52 after 7 shuffles lies in [0.334, 0.335), the published Bayer–Diaconis value | `native_decide` |
| `tv100_k7` … `tv100_k12` | **The Commander headline numbers**: exact TV distance at n = 100 lies in [0.747, 0.748), [0.428, 0.429), [0.223, 0.224), [0.112, 0.113), [0.056, 0.057), [0.028, 0.029) for k = 7…12 shuffles, computed in exact rational arithmetic from the formula and Eulerian counts | `native_decide` |

Because every (cut, interleaving) pair has probability exactly `a⁻ⁿ`
(Lemma 2 of `docs/PROOF.md`), equality of the *counts* in the first three
rows is equality of the *distributions*: this is the Bayer–Diaconis
formula, verified permutation-by-permutation.

## Scope — what this does and does not establish

- `native_decide` theorems rest on Lean's native-code reduction axiom
  (`Lean.ofReduceBool`); the file ends with `#print axioms` reports so
  this is visible, not hidden. The `decide` theorems are pure kernel
  checks depending on no axioms at all.
- Lean core (4.35) no longer ships `Nat.choose` / `Nat.factorial`, so the
  file defines its own `choose` (multiplicative formula, exact division
  at each step) and `fact`.
- Not claimed: a symbolic, all-n formal proof of the counting argument
  inside Lean. That argument is proved on paper in `docs/PROOF.md`;
  formalizing it for general n (stars-and-bars as a Lean bijection) is
  future work.
