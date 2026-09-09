# GMRES, sGMRES, RR, and sRR in MATLAB

Concise MATLAB implementations of four Krylov-subspace methods:

- `gmres.m` -- fixed-dimension GMRES with full Arnoldi;
- `sgmres.m` -- sketched GMRES with truncated Arnoldi and an SRDCT;
- `rr.m` -- Rayleigh--Ritz extraction from a full Arnoldi basis;
- `srr.m` -- sketched Rayleigh--Ritz extraction from a truncated basis.

This first milestone contains the methods, shared numerical helpers, and
correctness tests. Timing sweeps and crossover heatmaps are intentionally
deferred until the source has been reviewed.

## Requirements

- MATLAB R2026a
- Signal Processing Toolbox, for the unitary `dct` implementation and its
  `Type=2` and `Type=4` options

The code uses matrices or function handles as linear operators. Sparse
matrices work without conversion to dense storage.

## Quick start

Run MATLAB from this repository, or add its root directory to the path:

```matlab
addpath("/path/to/cs6220-sGMRES-sRR-matlab")

n = 200;
A = gallery("grcar",n) + 2*speye(n);
b = ones(n,1);
v1 = (1:n)';
d = 40;

[x,gmresInfo] = gmres(A,b,d);
[xs,sgmresInfo] = sgmres(A,b,d,Truncation=4,Seed=1);

[theta,U,rrInfo] = rr(A,v1,d,NumEigenpairs=5);
[thetas,Us,srrInfo] = srr(A,v1,d, ...
    Truncation=2,NumEigenpairs=5,Seed=1);
```

A matrix-free operator has the same interface:

```matlab
applyA = @(x) A*x;
x = gmres(applyA,b,d);
theta = rr(applyA,v1,d,NumEigenpairs=5);
```

## Public interfaces

```matlab
[x,info] = gmres(A,b,d,InitialGuess=x0)

[x,info] = sgmres(A,b,d, ...
    InitialGuess=x0,SketchSize=s,Truncation=k,Seed=seed)

[theta,U,info] = rr(A,v1,d,NumEigenpairs=q)

[theta,U,info] = srr(A,v1,d, ...
    SketchSize=s,Truncation=k,NumEigenpairs=q,Seed=seed)
```

GMRES performs exactly `d` Arnoldi steps unless happy breakdown occurs. It
does not stop early at a residual tolerance; that fixed-work behavior will
make the later timing comparison easier to interpret.

The defaults follow the parameter choices used for the associated numerical
examples:

| Method | Arnoldi orthogonalization | Sketch |
|---|---|---|
| GMRES | all previous vectors, two passes | none |
| sGMRES | latest `k=4` vectors, two passes | DCT-II, `s=min(2(d+1),n)` |
| RR | all previous vectors, two passes | none |
| sRR | latest `k=2` vectors, two passes | DCT-IV, `s=min(4d,n)` |

Each sketched call constructs one seeded linear map and applies that same map
to every related vector or matrix. The random draw uses a local `RandStream`,
so it does not alter MATLAB's global random state.

The `info` structures report the achieved Krylov dimension, breakdown status,
true residual diagnostics, and method-specific quantities such as sketch size
and the condition number of the reduced QR factor. RR and sRR select Ritz
values with largest real part.

## Arnoldi implementation

MATLAB does not provide a public `arnoldi` function returning the basis and
Hessenberg matrix. The visible helpers
`krylov.arnoldi` and `krylov.truncatedArnoldi` implement the recurrence
directly with matrix multiplication, inner products, and norms.

`gallery("krylov",A,v1,d)` is not used: it forms the raw vectors
`[v1,A*v1,...,A^(d-1)*v1]` without orthogonalizing them.

The remaining shared helpers are also visible in the `+krylov` package:

- `krylov.srdct` constructs and applies a signed, subsampled unitary DCT;
- `krylov.asOperator` gives matrices and function handles one interface;
- `krylov.pivotedQR` detects rank-deficient sketched reduced problems;
- `krylov.ritzPairs` sorts and reconstructs normalized Ritz pairs.

## A note about MATLAB's built-in `gmres`

This repository intentionally contains a function named `gmres.m`. While the
repository root is on the MATLAB path, it shadows MathWorks' function of the
same name. Use

```matlab
which gmres -all
```

to inspect path resolution. Our implementation is deliberate: MathWorks
`gmres` uses Householder transformations, Givens rotations, convergence and
stagnation logic, and other production-solver machinery. Using it as the
classical timing baseline would not isolate the algorithmic difference between
full and truncated Arnoldi.

## Tests

From the repository root:

```matlab
results = runtests("tests");
assertSuccess(results);
```

The tests check the Arnoldi relation, orthogonality, recent-vector truncated
orthogonalization, SRDCT isometry and reproducibility, matrix-free operators,
GMRES against a direct solve, RR against `eig`, and the limiting cases where a
full sketch and full truncation make the sketched methods agree with their
classical counterparts.

## AI acknowledgment

OpenAI Codex assisted with the initial implementation, tests, and
documentation in this repository. The author reviewed the resulting code and
is responsible for its correctness and presentation.

## Reference

The sketched methods and default parameters follow Yuji Nakatsukasa and Joel
A. Tropp, [Fast and Accurate Randomized Algorithms for Linear Systems and
Eigenvalue Problems](https://www.tropp.caltech.edu/papers/NT24-Fast-Accurate-SIMAX.pdf),
*SIAM Journal on Matrix Analysis and Applications*.
