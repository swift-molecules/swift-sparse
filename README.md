# swift-sparse

Use this package when a symmetric matrix is too large to store densely and most of its entries are zero: the stiffness and mass matrices a finite-element assembly produces, a graph Laplacian. Start from `Sparse.Symmetric<Scalar>`, built from `Entry` values the way element matrices are added into a global one. Multiply with `applied(to:)`, form `K − σM` with `adding(_:scaledBy:)`, and solve with `factorization()`, an `L D Lᵀ` factorization in a nested-dissection ordering whose `inertia` counts the matrix's negative eigenvalues. When many matrices share one pattern (each shift of one pencil), compute the `Sparse.Elimination` once and pass it to `Sparse.Symmetric.Factorization(_:elimination:)`. For the eigenvalues of a pencil `K x = λ M x`, use swift-eigen, which builds on this package. `Sparse Test Support` holds band, grid and scattered fixtures, the first two with closed-form spectra.

## Laws

- An entry adds its value to the entry at `{row, column}`, which is the same entry as `{column, row}`. Entries at the same place add up, and the place stays stored even when they cancel. The diagonal is always stored.
- Storage is compressed rows of the upper triangle: row `i` holds its columns `j ≥ i` in increasing order, the diagonal first. Two matrices are equal when their patterns and stored values are equal.
- `adding(_:scaledBy:)` stores the union of the two patterns, whatever the values, so every `K − σM` of one pencil has the same pattern and shares one `Sparse.Elimination`.
- `Sparse.Ordering.dissection(of:)` returns a permutation of every index. It numbers each piece's two sides before the level that cuts them, so eliminating one side never fills the other. A graph in several pieces is dissected piece by piece. On a 100 × 100 grid it leaves fewer than a third of the natural order's fill.
- The factorization takes each pivot in the ordering's place, without pivoting. A zero or non-finite pivot fails with `.singular`, naming the row in the matrix's own numbering. For a positive definite matrix every pivot is positive.
- `inertia` counts the signs of `D`, which by Sylvester's law are the signs of the matrix's eigenvalues. For `K − σM` with `M` positive definite, `inertia.negative` is the number of eigenvalues of `K x = λ M x` below `σ`.
- `solved(_:)` returns `x` with `A x = b`, up to the rounding of the factorization.

## Storage and ordering

Compressed rows rather than a skyline: a skyline stores every entry between the first nonzero of a row and its diagonal, so its cost is fixed by the ordering's profile. Compressed storage with an elimination tree stores only the fill the ordering produces, so the ordering can be changed without changing the storage. Nested dissection rather than a band (Cuthill–McKee) ordering: on an `n`-point plate or shell mesh a band ordering fills about `n^1.5` entries and a dissection about `n log n`. Every step of a shift-invert eigen solve is a pair of triangular solves with the factor, so that difference is paid thousands of times.
