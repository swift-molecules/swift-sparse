extension Sparse.Symmetric {
    /// `P A Pᵀ = L D Lᵀ` with `L` unit lower triangular and `D` diagonal, for the ordering `P` of an elimination.
    ///
    /// No pivoting: each pivot is taken in the ordering's place. That is exact for positive definite matrices and is
    /// how a shifted pencil `K − σM` is factored for its inertia; a zero or non-finite pivot fails with `.singular`.
    public struct Factorization {
        /// The elimination this factorization was computed on.
        public let elimination: Sparse.Elimination
        @usableFromInline
        internal let rows: [Int]
        @usableFromInline
        internal let values: [Scalar]
        @usableFromInline
        internal let pivots: [Scalar]
    }
}

extension Sparse.Symmetric.Factorization: Sendable where Scalar: Sendable {}

extension Sparse.Symmetric {
    /// The factorization of this matrix in a nested-dissection ordering.
    @inlinable
    public func factorization() throws(Sparse.Error) -> Factorization {
        try Factorization(self, elimination: Sparse.Elimination(pattern))
    }
}

extension Sparse.Symmetric.Factorization {
    /// Factors `matrix` on `elimination`, which must be an elimination of the matrix's pattern (`.pattern` otherwise).
    @inlinable
    public init(_ matrix: Sparse.Symmetric<Scalar>, elimination: Sparse.Elimination) throws(Sparse.Error) {
        guard matrix.pattern == elimination.pattern else { throw .pattern }
        let order = matrix.order
        let parent = elimination.parent
        let columnStarts = elimination.columnStarts
        let rowStarts = elimination.rowStarts
        let rowColumns = elimination.rowColumns
        let sources = elimination.sources
        let input = matrix.values
        var rows = [Int](repeating: 0, count: elimination.count)
        var values = [Scalar](repeating: 0, count: elimination.count)
        var pivots = [Scalar](repeating: 0, count: order)
        var filled = [Int](repeating: 0, count: order)
        var flag = [Int](repeating: -1, count: order)
        var path = [Int](repeating: 0, count: order)
        var y = [Scalar](repeating: 0, count: order)
        // Up-looking: row k of L solves a triangular system along the elimination tree (Davis's LDL).
        for k in 0..<order {
            var top = order
            flag[k] = k
            for stored in rowStarts[k]..<rowStarts[k + 1] {
                var i = rowColumns[stored]
                y[i] += input[sources[stored]]
                var length = 0
                while i < k && flag[i] != k {
                    path[length] = i
                    length += 1
                    flag[i] = k
                    i = parent[i]
                }
                while length > 0 {
                    length -= 1
                    top -= 1
                    path[top] = path[length]
                }
            }
            var pivot = y[k]
            y[k] = 0
            while top < order {
                let i = path[top]
                top += 1
                let yi = y[i]
                y[i] = 0
                let start = columnStarts[i]
                let end = start + filled[i]
                for p in start..<end {
                    y[rows[p]] -= values[p] * yi
                }
                let l = yi / pivots[i]
                pivot -= l * yi
                rows[end] = k
                values[end] = l
                filled[i] += 1
            }
            guard pivot != 0, pivot.isFinite else {
                throw .singular(elimination.ordering.permutation[k])
            }
            pivots[k] = pivot
        }
        self.elimination = elimination
        self.rows = rows
        self.values = values
        self.pivots = pivots
    }

    /// The order of the factored matrix.
    @inlinable
    public var order: Int { pivots.count }

    /// How many eigenvalues of the factored matrix are negative and how many positive: the signs of `D`.
    ///
    /// For `K − σM` with `M` positive definite, `negative` is the number of eigenvalues of `K x = λ M x` below `σ`.
    @inlinable
    public var inertia: Sparse.Inertia {
        var negative = 0
        for pivot in pivots where pivot < 0 { negative += 1 }
        return Sparse.Inertia(negative: negative, positive: pivots.count - negative)
    }

    /// The `x` with `A x = b`.
    @inlinable
    public func solved(_ b: [Scalar]) -> [Scalar] {
        var x = [Scalar](repeating: 0, count: order)
        solve(b, into: &x)
        return x
    }

    /// Writes the `x` with `A x = b` into `x`; both have the matrix's order.
    @inlinable
    public func solve(_ b: [Scalar], into x: inout [Scalar]) {
        precondition(b.count == order && x.count == order, "a vector of the matrix's order")
        let permutation = elimination.ordering.permutation
        let starts = elimination.columnStarts
        var z = [Scalar](repeating: 0, count: order)
        for k in 0..<order { z[k] = b[permutation[k]] }
        for column in 0..<order {
            let zc = z[column]
            if zc != 0 {
                for p in starts[column]..<starts[column + 1] {
                    z[rows[p]] -= values[p] * zc
                }
            }
        }
        for k in 0..<order { z[k] /= pivots[k] }
        var column = order - 1
        while column >= 0 {
            var sum = z[column]
            for p in starts[column]..<starts[column + 1] {
                sum -= values[p] * z[rows[p]]
            }
            z[column] = sum
            column -= 1
        }
        for k in 0..<order { x[permutation[k]] = z[k] }
    }
}
