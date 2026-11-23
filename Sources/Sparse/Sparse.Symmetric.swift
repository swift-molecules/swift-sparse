extension Sparse {
    /// A symmetric square matrix that stores only the entries of its upper triangle that were given.
    ///
    /// Storage is compressed rows of the upper triangle: for each row, the columns at or right of the
    /// diagonal and their values. The diagonal is always stored, as zero where no entry gave it.
    public struct Symmetric<Scalar: BinaryFloatingPoint> {
        /// Where the stored entries sit.
        public let pattern: Pattern
        @usableFromInline
        internal var values: [Scalar]

        @usableFromInline
        internal init(pattern: Pattern, values: [Scalar]) {
            self.pattern = pattern
            self.values = values
        }
    }
}

extension Sparse.Symmetric: Equatable {}
extension Sparse.Symmetric: Hashable where Scalar: Hashable {}
extension Sparse.Symmetric: Sendable where Scalar: Sendable {}

extension Sparse.Symmetric {
    /// One addition to a matrix: `value` added to the entry at `row`, `column`, which is also the entry at `column`, `row`.
    public struct Entry {
        public var row: Int
        public var column: Int
        public var value: Scalar

        @inlinable
        public init(_ row: Int, _ column: Int, _ value: Scalar) {
            self.row = row
            self.column = column
            self.value = value
        }
    }
}

extension Sparse.Symmetric.Entry: Equatable {}
extension Sparse.Symmetric.Entry: Hashable where Scalar: Hashable {}
extension Sparse.Symmetric.Entry: Sendable where Scalar: Sendable {}

extension Sparse.Symmetric {
    /// The `order × order` matrix that is the sum of `entries`, the way a finite-element assembly adds element
    /// matrices into a global one.
    ///
    /// Each entry adds its value to the entry at `{row, column}`, which is the same entry as `{column, row}`:
    /// give each off-diagonal pair once, from either triangle. Entries at the same place add up, and the place is
    /// stored even when they sum to zero. Fails with `.order` for a negative order and with `.index` for a row or
    /// column outside `0..<order`.
    @inlinable
    public init(order: Int, entries: some Sequence<Entry>) throws(Sparse.Error) {
        guard order >= 0 else { throw .order }
        var rows: [Int] = []
        var columns: [Int] = []
        var values: [Scalar] = []
        for entry in entries {
            guard entry.row >= 0, entry.row < order, entry.column >= 0, entry.column < order else {
                throw .index(row: entry.row, column: entry.column)
            }
            rows.append(Swift.min(entry.row, entry.column))
            columns.append(Swift.max(entry.row, entry.column))
            values.append(entry.value)
        }
        for index in 0..<order {
            rows.append(index)
            columns.append(index)
            values.append(0)
        }
        self = Self.compressed(order: order, rows: rows, columns: columns, values: values)
    }

    /// Sorts triplets by row then column with two counting passes, then adds up the duplicates.
    @usableFromInline
    internal static func compressed(order: Int, rows: [Int], columns: [Int], values: [Scalar]) -> Self {
        let total = rows.count
        var byColumn = [Int](repeating: 0, count: order + 1)
        for column in columns { byColumn[column + 1] += 1 }
        for index in 0..<order { byColumn[index + 1] += byColumn[index] }
        var orderByColumn = [Int](repeating: 0, count: total)
        for entry in 0..<total {
            orderByColumn[byColumn[columns[entry]]] = entry
            byColumn[columns[entry]] += 1
        }
        var byRow = [Int](repeating: 0, count: order + 1)
        for row in rows { byRow[row + 1] += 1 }
        for index in 0..<order { byRow[index + 1] += byRow[index] }
        var next = Array(byRow[0..<Swift.max(order, 0)])
        var sorted = [Int](repeating: 0, count: total)
        for entry in orderByColumn {
            sorted[next[rows[entry]]] = entry
            next[rows[entry]] += 1
        }
        var starts = [Int](repeating: 0, count: order + 1)
        var stored: [Int] = []
        var sums: [Scalar] = []
        stored.reserveCapacity(total)
        sums.reserveCapacity(total)
        for row in 0..<order {
            starts[row] = stored.count
            var last = -1
            for position in byRow[row]..<byRow[row + 1] {
                let entry = sorted[position]
                if columns[entry] == last {
                    sums[sums.count - 1] += values[entry]
                } else {
                    last = columns[entry]
                    stored.append(last)
                    sums.append(values[entry])
                }
            }
        }
        starts[order] = stored.count
        return Self(pattern: Sparse.Pattern(order: order, starts: starts, columns: stored), values: sums)
    }
}

extension Sparse.Symmetric {
    /// The number of rows, and of columns.
    @inlinable
    public var order: Int { pattern.order }

    /// The entry at `row`, `column`: zero where nothing is stored.
    @inlinable
    public subscript(row: Int, column: Int) -> Scalar {
        let low = Swift.min(row, column)
        let high = Swift.max(row, column)
        guard let position = pattern.position(row: low, column: high) else { return 0 }
        return values[position]
    }

    /// The stored entries of the upper triangle, row by row and left to right, the diagonal included.
    @inlinable
    public var entries: [Entry] {
        var result: [Entry] = []
        result.reserveCapacity(values.count)
        for row in 0..<order {
            for position in pattern.starts[row]..<pattern.starts[row + 1] {
                result.append(Entry(row, pattern.columns[position], values[position]))
            }
        }
        return result
    }

    /// The diagonal entries, first row first.
    @inlinable
    public var diagonal: [Scalar] {
        (0..<order).map { values[pattern.starts[$0]] }
    }

    /// The largest sum of absolute values along a row, which for a symmetric matrix is both its 1-norm and its ∞-norm.
    @inlinable
    public var norm: Scalar {
        var sums = [Scalar](repeating: 0, count: order)
        for row in 0..<order {
            for position in pattern.starts[row]..<pattern.starts[row + 1] {
                let size = Swift.abs(values[position])
                sums[row] += size
                let column = pattern.columns[position]
                if column != row { sums[column] += size }
            }
        }
        return sums.max() ?? 0
    }
}
