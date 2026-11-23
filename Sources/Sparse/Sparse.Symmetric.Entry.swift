extension Sparse.Symmetric {
    /// One addition to a matrix: `value` added to the entry at `row`, `column`, which is also the entry at `column`, `row`.
    public struct Entry {
        public var row: Ordinal
        public var column: Ordinal
        public var value: Scalar

        @inlinable
        public init(_ row: Ordinal, _ column: Ordinal, _ value: Scalar) {
            self.row = row
            self.column = column
            self.value = value
        }
    }
}

extension Sparse.Symmetric.Entry: Equatable {}
extension Sparse.Symmetric.Entry: Hashable where Scalar: Hashable {}
extension Sparse.Symmetric.Entry: Sendable where Scalar: Sendable {}
