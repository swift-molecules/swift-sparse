extension Sparse {
    /// The symbolic half of a factorization: an ordering, the elimination tree it gives, and where every entry of the
    /// factor will sit. It depends only on a pattern, so one elimination serves every matrix on that pattern, such as
    /// `K − σM` for each shift `σ`.
    public struct Elimination: Hashable, Sendable {
        /// The pattern this elimination factors.
        public let pattern: Pattern
        /// The ordering the factor is computed in.
        public let ordering: Ordering
        /// `parent[k]` is the parent of position `k` in the elimination tree, or -1 at a root.
        @usableFromInline
        internal let parent: [Int]
        /// Where each column of the factor starts in its storage; column `k` holds the rows below position `k`.
        @usableFromInline
        internal let columnStarts: [Int]
        /// The stored entries of the permuted lower triangle, row by row: their columns and where their values come from.
        @usableFromInline
        internal let rowStarts: [Int]
        @usableFromInline
        internal let rowColumns: [Int]
        @usableFromInline
        internal let sources: [Int]
    }
}

extension Sparse.Elimination {
    /// The elimination of `pattern` in a nested-dissection ordering.
    public init(_ pattern: Sparse.Pattern) {
        self.init(pattern, ordering: .dissection(of: pattern))
    }

    /// The elimination of `pattern` in `ordering`, which must have the pattern's order.
    public init(_ pattern: Sparse.Pattern, ordering: Sparse.Ordering) {
        let order = pattern.order
        precondition(ordering.permutation.count == order, "an ordering of the pattern's order")
        let position = ordering.position
        // WHY: the permuted lower triangle, row by row, is what the up-looking factorization reads.
        var rowStarts = [Int](repeating: 0, count: order + 1)
        var row = 0
        while row < order {
            var stored = pattern.starts[row]
            let storedEnd = pattern.starts[row + 1]
            while stored < storedEnd {
                let column = pattern.columns[stored]
                rowStarts[Swift.max(position[row], position[column]) + 1] += 1
                stored += 1
            }
            row += 1
        }
        for index in 0..<order { rowStarts[index + 1] += rowStarts[index] }
        var next = Array(rowStarts[0..<order])
        var rowColumns = [Int](repeating: 0, count: pattern.count)
        var sources = [Int](repeating: 0, count: pattern.count)
        row = 0
        while row < order {
            var stored = pattern.starts[row]
            let storedEnd = pattern.starts[row + 1]
            while stored < storedEnd {
                let column = pattern.columns[stored]
                let a = position[row]
                let b = position[column]
                let high = Swift.max(a, b)
                rowColumns[next[high]] = Swift.min(a, b)
                sources[next[high]] = stored
                next[high] += 1
                stored += 1
            }
            row += 1
        }
        // WHY: the elimination tree and each column count of the factor follow Liu, as in Davis's LDL.
        var parent = [Int](repeating: -1, count: order)
        var flag = [Int](repeating: -1, count: order)
        var counts = [Int](repeating: 0, count: order)
        var k = 0
        while k < order {
            flag[k] = k
            var stored = rowStarts[k]
            let storedEnd = rowStarts[k + 1]
            while stored < storedEnd {
                var i = rowColumns[stored]
                while i < k && flag[i] != k {
                    if parent[i] == -1 { parent[i] = k }
                    counts[i] += 1
                    flag[i] = k
                    i = parent[i]
                }
                stored += 1
            }
            k += 1
        }
        var columnStarts = [Int](repeating: 0, count: order + 1)
        for k in 0..<order { columnStarts[k + 1] = columnStarts[k] + counts[k] }
        self.pattern = pattern
        self.ordering = ordering
        self.parent = parent
        self.columnStarts = columnStarts
        self.rowStarts = rowStarts
        self.rowColumns = rowColumns
        self.sources = sources
    }

    /// The number of entries of the factor below its diagonal: the fill the ordering leaves.
    @inlinable
    public var count: Int { columnStarts.last ?? 0 }
}
