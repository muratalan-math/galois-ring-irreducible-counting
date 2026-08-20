# Reproducibility code for irreducible-polynomial counting over Galois rings

This repository contains the computational verification code accompanying the manuscript

> **Effective Enumeration of Irreducible Polynomials over Galois Rings: Low-Degree Lifting Automata and Characteristic \(p^3\)**

The manuscript develops explicit counting formulas and finite lifting automata for monic irreducible polynomials over Galois rings. **No computer calculation is used as part of a mathematical proof.** The code in this repository is a reproducibility and verification artifact: it independently checks selected finite cases, finite polynomial maps/fibers, theorem boundary cases, and symbolic recurrence identities.

## Authors

- **Murat Alan** (corresponding author), Department of Mathematics, Yildiz Technical University; ORCID 0000-0003-2031-2725.
- **Kadriye Şimşek Alan**, Department of Mathematical Engineering, Yildiz Technical University; ORCID 0000-0001-6751-8013.

## Quick start

### Requirements

- a C++17 compiler (reference run: `g++ 14.2.0`)
- Python 3 (reference run: `Python 3.13.5`)
- SymPy 1.14.0

Install the Python dependency with:

```bash
python3 -m pip install -r requirements.txt
```

Run the complete verification suite with:

```bash
bash run_all.sh
```

A successful reference run ends with:

```text
ALL REPRODUCIBILITY CHECKS PASSED
```

On the archived reference environment, the full suite takes roughly tens of seconds; the weighted-origin enumeration is the heaviest individual component. Runtime is not used as evidence for any theorem or complexity claim.

## What is checked

The repository deliberately separates three logically different kinds of verification.

### 1. `ground_truth/` — direct finite-ring enumeration

These programs **do not call the manuscript recurrences**.

- `z8_local_all_degrees.cpp` directly enumerates all local monic polynomials over `Z/8Z` reducing to `x^e` for `1 <= e <= 6` and marks all nontrivial monic factor products; this independently reproduces the full characteristic-`p^3` benchmark `J_{2,3}(1),...,J_{2,3}(6)`.
- `local_direct_enumeration.cpp` constructs the relevant finite rings explicitly, including non-prime residue-field examples, enumerates root-bearing polynomials and monic factor products, and counts the complement. For quartics it additionally reports the root-free count \(U_s\) and the root-free-with-quadratic-factor count \(B_s=Q^{s-1}D_s\) separately, so both automaton outputs are checked against ground truth (see the \(s=4\), \(Q=5\) point).
- `global_z8_direct.cpp` directly enumerates monic irreducibles over \(\mathbb Z/8\mathbb Z\) in degrees 1–4, independently of the structural decomposition formula.

### 2. `finite_maps/` and `proof_boundaries/` — direct finite-map/fiber checks

These programs enumerate the finite polynomial maps occurring in the quartic factor automaton and test delicate boundary mechanisms, including mixed-precision fibers, the nonzero repeated-quadratic state, weighted origin contributions, the \(\mathbb F_{25}\) residual classification, valuation inequalities, and the exact \(Q^2\) overlap fiber at the boundary.

### 3. `symbolic/` — formula and recurrence consistency

These scripts intentionally re-evaluate formulas appearing in the manuscript. They are **not** advertised as independent ground truth.

- `denominator_check.py` expands the announced factor denominator and verifies every coefficient.
- `factor_start_index_check.py` reconstructs the factor sequence from the exact decomposition (without using the order-12 recurrence), verifies the nonzero \(T^{12}\) numerator coefficient, and checks vanishing convolutions from level 13 onward.
- `closed_formula_checks.py` evaluates selected closed formulas and recurrence values.
- `quartic_automaton.py` assembles the full tame-quartic algorithm of Section 10 in **two algebraically equivalent implementations** — the scalar finite-state recurrences and the fixed \(9\times9\) root / \(12\times12\) factor transition matrices via binary powering — checks the two paths agree for \(Q\in\{5,7,25\}\) and \(s\le40\), and cross-checks the assembled \((U_s,D_s,J)\) against the direct ground-truth values (including the \(s=4\) point) and the characteristic-\(p^3\) closed form.

See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for a theorem-to-code map.

## Repository layout

```text
ground_truth/       direct finite-ring enumeration
finite_maps/        direct finite polynomial-map and fiber checks
proof_boundaries/   finite checks of delicate proof boundary cases
symbolic/           symbolic formula/recurrence consistency checks
expected/           archived deterministic reference outputs
.github/workflows/  GitHub Actions reproducibility workflow
```

Temporary build products are written to `_build/` and `_run_outputs/`; both are ignored by Git. The archived repository release version is recorded in `VERSION` (current submission artifact: `1.2.1`).

## Expected outputs and integrity

Every deterministic test output is archived under `expected/`. The runner compares fresh outputs with these files and stops on any discrepancy.

For an exact release, repository-file hashes are recorded in `SHA256SUMS.txt`. Verify them with:

```bash
bash verify_checksums.sh
```

The checksum file intentionally excludes itself.

## Scope and interpretation

The quartic automaton tests are restricted to the proved tame range: the end-to-end driver is run only at `Q=5,7,25`, corresponding to residue characteristics at least five. The `Q=2` and `Q=3` checks elsewhere in the repository concern the degree-two, cubic, and characteristic-`p^3` results, which are valid in those characteristics. **No excluded-characteristic quartic continuation experiment is part of this release.**

Direct enumeration is restricted to finite parameter ranges where exhaustive computation is practical. Agreement with those computations is strong evidence against normalization, indexing, and fiber-multiplicity errors, but it is not substituted for proof. The manuscript's effective evaluation and bit-complexity results come from proved finite-state recurrences and fixed-size transition matrices.

## Reproducibility policy

- Ground-truth enumeration is kept separate from formula/recurrence checking.
- Expected outputs are plain text and human-readable.
- Tool versions are recorded in `TOOL_VERSIONS.txt`.
- The complete suite is run automatically by GitHub Actions.
- Release hashes identify the exact archived artifact.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). When the manuscript receives a DOI or arXiv identifier, the citation metadata can be updated without changing the verification code.

## License

The source code in this repository is released under the MIT License; see [`LICENSE`](LICENSE).
