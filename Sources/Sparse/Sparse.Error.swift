extension Sparse {
    /// Why a sparse matrix could not be built, combined or factored.
    public enum Error: Swift.Error, Hashable, Sendable {
        /// An order is negative, or two matrices that must have the same order do not.
        case order
        /// An entry names a row or column outside the matrix.
        case index(row: Int, column: Int)
        /// A matrix has a different pattern from the elimination it is factored with.
        case pattern
        /// The pivot of the row at this index, in the matrix's own numbering, is zero or not finite.
        case singular(Int)
    }
}
