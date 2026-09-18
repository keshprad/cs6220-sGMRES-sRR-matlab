# GMRES, sGMRES, RR, and sRR in MATLAB

Concise MATLAB implementations of four Krylov-subspace methods:

- `gmres.m` -- fixed-dimension GMRES with full Arnoldi;
- `sgmres.m` -- sketched GMRES with truncated Arnoldi and an SRDCT;
- `rr.m` -- Rayleigh--Ritz extraction from a full Arnoldi basis;
- `srr.m` -- sketched Rayleigh--Ritz extraction from a truncated basis.

The repository also contains a MATLAB-native paired timing experiment that
produces raw data, summaries, and crossover line charts.

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
k = 40;
ellGmres = 4;
ellRr = 2;

[x,gmresInfo] = gmres(A,b,k);
[xs,sgmresInfo] = sgmres(A,b,k,Truncation=ellGmres,Seed=1);

[theta,U,rrInfo] = rr(A,v1,k,NumEigenpairs=5);
[thetas,Us,srrInfo] = srr(A,v1,k, ...
    Truncation=ellRr,NumEigenpairs=5,Seed=1);
```

A matrix-free operator has the same interface:

```matlab
applyA = @(x) A*x;
x = gmres(applyA,b,k);
theta = rr(applyA,v1,k,NumEigenpairs=5);
```

## Crossover experiment

Run a short check or the normal lecture-development sweep from the repository
root:

```matlab
[raw,summary,paths] = experiments.runCrossover( ...
    Preset="smoke",OutputDirectory="results/smoke");

[raw,summary,paths] = experiments.runCrossover( ...
    Preset="quick",OutputDirectory="results/quick");
```

The larger `Preset="lecture"` sweep uses grid sizes 64 and 128, Krylov
dimensions 16 through 256, six paired repetitions, and five independent
sketch trials. Existing nonempty output directories are protected. Pass
`Force=true` only when deliberately replacing the named result bundle.

GMRES and sGMRES use the sparse five-point matrix returned by
`gallery("poisson",gridSize)`. A seeded local random stream constructs a
normalized exact solution and the right-hand side is `b=A*xExact`; MATLAB's
global random state is unchanged. RR and sRR use the same explicit sparse
Poisson matrix, with the seeded normalized vector as the Arnoldi start.

Every bundle contains `raw.csv`, `summary.csv`, and `metadata.json`, plus
`gmres_crossover.png/.pdf` and `rr_crossover.png/.pdf` when plotting is
enabled. Each chart has one line per grid size, with Krylov dimension on a
logarithmic horizontal axis. The left panel plots the median paired
`classical core time / sketched core time` on a logarithmic vertical axis;
values above one favor the sketched method. Error bars show the 25th–75th
percentiles of the same paired ratios used for the median (30 pairs per
setting in the lecture sweep). They show timing spread, not confidence
intervals. The right panel plots the median
sketched/classical true-residual ratio on a linear vertical axis, with:

- a shaded band for the 25th–75th percentiles;
- whiskers spanning the minimum and maximum;
- faint points for the independent sketch trials.

Residual spread uses one median paired residual ratio per sketch trial, so
repeated timings do not count as independent residual samples. Quartiles use
linear interpolation at ranks `1+(N-1)*p`. The bands describe observed spread,
not confidence intervals. Dashed lines mark a ratio of one. All plotted
settings use filled markers. Settings that fail the fixed-work or
accuracy gates are left as gaps. Explanatory text belongs in the figure
caption; it is not embedded beneath the chart.

To regenerate charts from saved measurements without rerunning the experiment:

```matlab
folder = "results/lecture";
raw = readtable(fullfile(folder,"raw.csv"),TextType="string");
metadata = jsondecode(fileread(fullfile(folder,"metadata.json")));
summary = experiments.summarize(raw,metadata.Config);
[paths,plotData] = experiments.plotCrossover(raw,summary,folder);
```

Recomputing the summary adds timing quartiles to older result bundles without
rerunning any timed method. `plotData` exposes the plotted medians, quartiles, extrema, and individual
sketch ratios for verification. The plotting interface now requires `raw`
before `summary`; old summary-only calls must be updated because a median
cannot reconstruct the spread.

The calls are paired and their order alternates to reduce timing bias.
Problem construction, file output, plotting, and the built-in `eigs` oracle
used to check RR targets are outside the method timings. The sketch draw and
DCT applications remain inside the sketched timings.

## Public interfaces

```matlab
[x,info] = gmres(A,b,k,InitialGuess=x0)

[x,info] = sgmres(A,b,k, ...
    InitialGuess=x0,SketchSize=s,Truncation=ell,Seed=seed)

[theta,U,info] = rr(A,v1,k,NumEigenpairs=q)

[theta,U,info] = srr(A,v1,k, ...
    SketchSize=s,Truncation=ell,NumEigenpairs=q,Seed=seed)
```

The notation matches the lecture notes: `k` is the Krylov dimension and
effective iteration count, while `ell` is the number of recent basis vectors
used by truncated Arnoldi.

GMRES performs exactly `k` Arnoldi steps unless happy breakdown occurs. It
does not stop early at a residual tolerance; that fixed-work behavior will
make the later timing comparison easier to interpret.

The defaults follow the parameter choices used for the associated numerical
examples:

| Method | Arnoldi orthogonalization | Sketch |
|---|---|---|
| GMRES | all previous vectors, two passes | none |
| sGMRES | latest `ell=4` vectors, two passes | DCT-II, `s=min(2(k+1),n)` |
| RR | all previous vectors, two passes | none |
| sRR | latest `ell=2` vectors, two passes | DCT-IV, `s=min(4k,n)` |

Each sketched call constructs one seeded linear map and applies that same map
to every related vector or matrix. The random draw uses a local `RandStream`,
so it does not alter MATLAB's global random state.

The `info` structures report the achieved Krylov dimension, breakdown status,
true residual diagnostics, `CoreTime`, `WallTime`, and method-specific
quantities such as sketch size and the condition number of the reduced QR
factor. RR and sRR select Ritz values with largest real part.

## Arnoldi implementation

MATLAB does not provide a public `arnoldi` function returning the basis and
Hessenberg matrix. The visible helpers
`krylov.arnoldi` and `krylov.truncatedArnoldi` implement the recurrence
directly with matrix multiplication, inner products, and norms.

`gallery("krylov",A,v1,k)` is not used: it forms the raw vectors
`[v1,A*v1,...,A^(k-1)*v1]` without orthogonalizing them.

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
