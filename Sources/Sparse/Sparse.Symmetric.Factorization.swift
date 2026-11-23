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
        var rows = [Int](repeating: 0, count: elimination.count)
        var values = [Scalar](repeating: 0, count: elimination.count)
        var pivots = [Scalar](repeating: 0, count: order)
        let failed = Self.factor(
            matrix.values,
            elimination: elimination,
            rows: &rows,
            values: &values,
            pivots: &pivots
        )
        if failed >= 0 { throw .singular(elimination.ordering.permutation[failed]) }
        self.elimination = elimination
        self.rows = rows
        self.values = values
        self.pivots = pivots
    }

    /// The numeric factorization, up-looking along the elimination tree (Davis's LDL): row `k` of `L` solves a sparse
    /// triangular system whose pattern is the tree path from each entry of row `k`. Returns the position of the first
    /// zero or non-finite pivot, or −1.
    @inlinable
    static func factor(
        _ input: [Scalar],
        elimination: Sparse.Elimination,
        rows: inout [Int],
        values: inout [Scalar],
        pivots: inout [Scalar]
    ) -> Int {
        let order = pivots.count
        let filled = UnsafeMutablePointer<Int>.allocate(capacity: Swift.max(order, 1))
        let flag = UnsafeMutablePointer<Int>.allocate(capacity: Swift.max(order, 1))
        let path = UnsafeMutablePointer<Int>.allocate(capacity: Swift.max(order, 1))
        let y = UnsafeMutablePointer<Scalar>.allocate(capacity: Swift.max(order, 1))
        defer {
            unsafe filled.deallocate()
            unsafe flag.deallocate()
            unsafe path.deallocate()
            unsafe y.deallocate()
        }
        unsafe filled.initialize(repeating: 0, count: order)
        unsafe flag.initialize(repeating: -1, count: order)
        unsafe y.initialize(repeating: 0, count: order)
        return input.withUnsafeBufferPointer { input in
            elimination.parent.withUnsafeBufferPointer { parent in
                elimination.columnStarts.withUnsafeBufferPointer { columnStarts in
                    elimination.rowStarts.withUnsafeBufferPointer { rowStarts in
                        elimination.rowColumns.withUnsafeBufferPointer { rowColumns in
                            elimination.sources.withUnsafeBufferPointer { sources in
                                rows.withUnsafeMutableBufferPointer { rows in
                                    values.withUnsafeMutableBufferPointer { values in
                                        pivots.withUnsafeMutableBufferPointer { pivots in
                                            var k = 0
                                            while k < order {
                                                var top = order
                                                unsafe flag[k] = k
                                                var stored = unsafe rowStarts[k]
                                                let storedEnd = unsafe rowStarts[k + 1]
                                                while stored < storedEnd {
                                                    var i = unsafe rowColumns[stored]
                                                    unsafe y[i] += input[sources[stored]]
                                                    var length = 0
                                                    while unsafe i < k && flag[i] != k {
                                                        unsafe path[length] = i
                                                        length += 1
                                                        unsafe flag[i] = k
                                                        i = unsafe parent[i]
                                                    }
                                                    while length > 0 {
                                                        length -= 1
                                                        top -= 1
                                                        unsafe path[top] = path[length]
                                                    }
                                                    stored += 1
                                                }
                                                var pivot = unsafe y[k]
                                                unsafe y[k] = 0
                                                while top < order {
                                                    let i = unsafe path[top]
                                                    top += 1
                                                    let yi = unsafe y[i]
                                                    unsafe y[i] = 0
                                                    let start = unsafe columnStarts[i]
                                                    let end = unsafe start + filled[i]
                                                    var p = start
                                                    while p < end {
                                                        unsafe y[rows[p]] -= values[p] * yi
                                                        p += 1
                                                    }
                                                    let l = unsafe yi / pivots[i]
                                                    pivot -= l * yi
                                                    unsafe rows[end] = k
                                                    unsafe values[end] = l
                                                    unsafe filled[i] += 1
                                                }
                                                guard pivot != 0, pivot.isFinite else { return k }
                                                unsafe pivots[k] = pivot
                                                k += 1
                                            }
                                            return -1
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
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
        let order = self.order
        let z = UnsafeMutablePointer<Scalar>.allocate(capacity: Swift.max(order, 1))
        defer { unsafe z.deallocate() }
        elimination.ordering.permutation.withUnsafeBufferPointer { permutation in
            elimination.columnStarts.withUnsafeBufferPointer { starts in
                rows.withUnsafeBufferPointer { rows in
                    values.withUnsafeBufferPointer { values in
                        pivots.withUnsafeBufferPointer { pivots in
                            b.withUnsafeBufferPointer { b in
                                var k = 0
                                while k < order {
                                    unsafe (z + k).initialize(to: b[permutation[k]])
                                    k += 1
                                }
                                var column = 0
                                while column < order {
                                    let zc = unsafe z[column]
                                    if zc != 0 {
                                        var p = unsafe starts[column]
                                        let end = unsafe starts[column + 1]
                                        while p < end {
                                            unsafe z[rows[p]] -= values[p] * zc
                                            p += 1
                                        }
                                    }
                                    column += 1
                                }
                                k = 0
                                while k < order {
                                    unsafe z[k] /= pivots[k]
                                    k += 1
                                }
                                column = order - 1
                                while column >= 0 {
                                    var sum = unsafe z[column]
                                    var p = unsafe starts[column]
                                    let end = unsafe starts[column + 1]
                                    while p < end {
                                        unsafe sum -= values[p] * z[rows[p]]
                                        p += 1
                                    }
                                    unsafe z[column] = sum
                                    column -= 1
                                }
                            }
                        }
                    }
                }
            }
        }
        elimination.ordering.permutation.withUnsafeBufferPointer { permutation in
            var k = 0
            while k < order {
                unsafe x[permutation[k]] = z[k]
                k += 1
            }
        }
    }
}
