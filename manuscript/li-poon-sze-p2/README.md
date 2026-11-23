# Elementary-symmetric realignment separation

For every unequal pair of local dimensions at least two, the note proves
`B_l < Btilde_l` simultaneously at every degree from two through the square
of the smaller dimension. It extends Li–Poon–Sze's equal-dimensional answer
(2011, Theorem 2.4) to their open problem (P2), arXiv:1103.2499v2, p. 4.
The original question is Lupo–Aniello–Scardicchio (2008), Section 6, p. 15,
item (ii).

The explicit witnesses also give positive-definite white-noise mixtures
strictly inside the CCNR constraint, with a positive mixing interval stated
in the note. Their normalized Choi matrices also give explicit non-entanglement-
breaking channels satisfying the first-degree singular-value bound, yet
violating every higher-degree separable profile bound.
Equal-dimensional equality, the sharp dimension threshold for
the naive CCNR bound, the two-qubit values, and the known exact CCNR ranges
are credited to their original sources. Exact remaining maxima are separated
from the strict-comparison result.

## Files and rebuilding

- `note.tex`: complete written argument, coverage table and bibliography.
- `note.pdf`: the mathematical note.
- `verify.py`, `verification.json`: the exact checks and their expected output.
- `metadata.json`: scientific scope, creators, license and pinned formal sources.
- `manifest.json`: SHA256 and byte counts of the components.

Build with Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, run `pdflatex note.tex` twice. The bibliography is in the
TeX source; there are no external figures or BibTeX dependencies.

## Exact checks

Use Python 3 with SymPy 1.14.0:

```sh
python3 verify.py > verification.json
```

The arithmetic is symbolic and rational, with outward square-root enclosures
at 70 decimal places. The script checks twelve complete realignment row-Gram
examples: six block states and their six white-noise mixtures at weight
`1/1000`. It also checks trace and noise identities, profile derivatives for
sizes 4–17 at every degree, 42,826 dimension/sector pairs with `2 <= m <= 50`,
the maximally mixed Alice marginals, and the exact two-qubit values.
These finite examples illustrate the identities;
the universal conclusions have written proofs in the note.

## Formal scope

The pinned
[`LiPoonSzeP2.result`](https://github.com/the-omega-institute/trureturing/blob/d13f56f01c2e4710683d27a74c67cba9ec344e19/D5/S3/Quantum/Entanglement/Witness/LiPoonSzeP2.lean#L806)
has quantifiers `forall m,n,l exists rho`: at every allowed degree it specifies
a complex density matrix satisfying CCNR and exceeding an upper bound valid
for every separable state. In that file, `main_witness`, `supported_witness`
and `uniform_witness` specify the spectra independently of the degree; the
written proof uses them to give one simultaneous witness.

[`RealignmentSchurBound.separable_esymm_bound`](https://github.com/the-omega-institute/trureturing/blob/d13f56f01c2e4710683d27a74c67cba9ec344e19/D5/S3/Quantum/Entanglement/Witness/RealignmentSchurBound.lean#L924)
gives the separable elementary-symmetric profile bound.
[`CorrelationMinorNormCovariance.actual_covariance_nuclear`](https://github.com/the-omega-institute/trureturing/blob/d13f56f01c2e4710683d27a74c67cba9ec344e19/D5/S3/Quantum/Entanglement/Witness/CorrelationMinorNormCovariance.lean#L162)
and `correlation_tail_purity` give the covariance and marginal-purity tail
inequalities in orthonormal Hermitian bases. The positive-definite white-noise
extension, channel consequence and equal-dimensional averaging argument are written proofs.
`metadata.json` identifies the exact source bytes and declarations.

## Creators and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang.
Creative Commons Attribution 4.0 International (CC BY 4.0).
