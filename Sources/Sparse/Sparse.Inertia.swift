extension Sparse {
    /// How many eigenvalues of a symmetric matrix are negative and how many are positive (Sylvester's law of inertia).
    public struct Inertia: Hashable, Sendable {
        /// The number of negative eigenvalues.
        public let negative: Int
        /// The number of positive eigenvalues.
        public let positive: Int

        @inlinable
        public init(negative: Int, positive: Int) {
            self.negative = negative
            self.positive = positive
        }
    }
}
