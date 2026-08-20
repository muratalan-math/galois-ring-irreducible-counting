# Reproducibility map

Accompanying manuscript by **Murat Alan** and **Kadriye Şimşek Alan**.

This document maps the principal computational checks to the corresponding mathematical claims in the accompanying manuscript. The checks are verification material only; none is used as a premise in a proof.

| Manuscript result / mechanism | Verification artifact | Verification type |
|---|---|---|
| Structural global-to-local decomposition (Theorem 3.1), small global cases | `ground_truth/global_z8_direct.cpp` | direct global finite-ring enumeration |
| Full `Z/8Z` characteristic-`p^3` benchmark `J_{2,3}(1),...,J_{2,3}(6)` | `ground_truth/z8_local_all_degrees.cpp` | direct enumeration of all local monic polynomials and all nontrivial monic factor products |
| Degree-two and characteristic-\(p^3\) local counts, selected parameters | `ground_truth/local_direct_enumeration.cpp` | direct local finite-ring enumeration |
| Cubic all-level recurrence, selected small levels in `Q=2,3` | `ground_truth/local_direct_enumeration.cpp`, `expected/cubic_direct.txt` | direct local enumeration at independent levels |
| Non-prime residue-field checks (e.g. \(\operatorname{GR}(2^s,2)\), \(\operatorname{GR}(5^2,2)\)) | `ground_truth/local_direct_enumeration.cpp` | explicit finite-ring arithmetic |
| Quartic mixed-precision coprime fiber (Lemma 8.14) | `finite_maps/mixed_precision_fibers.cpp` | direct finite-map/fiber enumeration |
| Nonzero repeated-quadratic state and refined digits (Lemma 8.8) | `finite_maps/nonzero_curve_state.cpp` | direct finite-image enumeration |
| Five residual origin types over \(\mathbb F_{25}\) | `finite_maps/gf25_residual_classification.py` | direct residual-product classification |
| Weighted origin contributions / corrected origin formula (Theorem 8.29), selected low levels | `finite_maps/weighted_origin.cpp` | direct weighted finite-image enumeration |
| Valuation inequalities and \(s=3\) regular-overlap boundary fiber | `proof_boundaries/proof_boundary_checks.py` | finite/exhaustive boundary check |
| Factor denominator and listed order-12 coefficients (Theorem 8.36) | `symbolic/denominator_check.py` | exact symbolic identity |
| Exact factor-recurrence start index \(s=13\) | `symbolic/factor_start_index_check.py` | exact symbolic reconstruction from factor decomposition |
| Selected degree-two, characteristic-\(p^3\), and cubic closed-form values | `symbolic/closed_formula_checks.py` | formula/recurrence consistency |
| End-to-end tame quartic count \(J_{Q,s}(4)\): scalar recurrences vs. fixed \(9\times9\)/\(12\times12\) matrices (Thm 10.2, Cor. 8.37) | `symbolic/quartic_automaton.py` | scalar-path vs. matrix-path consistency + cross-check against ground truth |
| Integrated quartic ground truth at \(s=4\): \(U_4,D_4,J_{5,4}(4)\) enumerated separately | `ground_truth/local_direct_enumeration.cpp` (args `5 4 1 4`) | direct finite-ring enumeration of root-free / quadratic-factor / irreducible counts |

## Independence hierarchy

The strongest numerical checks are those in `ground_truth/`: they enumerate finite rings directly and do not call the recurrences whose values they reproduce.

The `finite_maps/` programs sit one level closer to the proof: they enumerate the concrete finite polynomial maps and fibers used to derive the quartic automaton.

The `symbolic/` scripts are intentionally different: they confirm algebraic consistency of displayed formulas and are **not** called independent ground truth.

## Tame-quartic scope

No `p=2` or `p=3` quartic continuation experiment is included in this submission artifact. Quartic end-to-end checks use only `Q=5,7,25`; appearances of `Q=2` or `Q=3` belong to the all-characteristic degree-two/cubic/characteristic-`p^3` results.
