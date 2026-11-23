extension Sparse {
    /// A renumbering of a matrix's rows and columns, chosen so that its factor fills in little.
    public struct Ordering: Hashable, Sendable {
        /// `permutation[k]` is the original index placed at position `k`.
        public let permutation: [Int]
        /// `position[i]` is the position of original index `i`: the inverse of `permutation`.
        public let position: [Int]

        @usableFromInline
        internal init(permutation: [Int]) {
            var position = [Int](repeating: 0, count: permutation.count)
            permutation.enumerated().forEach { index, original in position[original] = index }
            self.permutation = permutation
            self.position = position
        }
    }
}

extension Sparse.Ordering {
    /// The identity: every index keeps its place.
    @inlinable
    public static func natural(order: Cardinal) -> Self {
        Self(permutation: Array(0..<Int(bitPattern: order)))
    }

    /// A nested dissection of the pattern's graph.
    ///
    /// Each connected piece larger than `leaf` vertices is cut by one level of a breadth-first search from a
    /// pseudo-peripheral vertex, the smallest level near the middle, trimmed of vertices that touch only one side.
    /// The two sides are numbered first, each dissected the same way, and the cut last, so eliminating a side never
    /// fills the other. Pieces of at most `leaf` vertices are numbered in reverse Cuthill–McKee order.
    public static func dissection(of pattern: Sparse.Pattern, leaf: Cardinal = 48) -> Self {
        var dissection = Dissection(pattern: pattern, leaf: Swift.max(Int(bitPattern: leaf), 1))
        return dissection.ordering()
    }
}
