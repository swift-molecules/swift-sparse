extension Sparse {
    /// Where the stored entries of a symmetric matrix sit: the upper triangle, row by row.
    ///
    /// Row `i` holds the columns `j >= i` of its stored entries in increasing order, and always holds
    /// its diagonal `i` first.
    public struct Pattern: Hashable, Sendable {
        /// The number of rows, and of columns.
        public let order: Int
        @usableFromInline
        internal let starts: [Int]
        @usableFromInline
        internal let columns: [Int]

        @usableFromInline
        internal init(order: Int, starts: [Int], columns: [Int]) {
            self.order = order
            self.starts = starts
            self.columns = columns
        }
    }
}

extension Sparse.Pattern {
    /// The number of stored entries in the upper triangle, the diagonal included.
    @inlinable
    public var count: Int { columns.count }

    /// The stored columns of row `row`, each at or right of the diagonal, in increasing order.
    @inlinable
    public func columns(in row: Int) -> ArraySlice<Int> {
        columns[starts[row]..<starts[row + 1]]
    }

    /// The position of the entry at `row`, `column` (with `row <= column`) among the stored entries, if stored.
    @inlinable
    internal func position(row: Int, column: Int) -> Int? {
        var low = starts[row]
        var high = starts[row + 1]
        while low < high {
            let middle = (low + high) / 2
            let found = columns[middle]
            if found == column { return middle }
            if found < column { low = middle + 1 } else { high = middle }
        }
        return nil
    }

    /// Every stored neighbour of every row, both triangles, the diagonal left out.
    @usableFromInline
    internal func adjacency() -> (starts: [Int], neighbours: [Int]) {
        var degree = [Int](repeating: 0, count: order + 1)
        var row = 0
        while row < order {
            var position = starts[row] + 1
            let positionEnd = starts[row + 1]
            while position < positionEnd {
                degree[row + 1] += 1
                degree[columns[position] + 1] += 1
                position += 1
            }
            row += 1
        }
        for row in 0..<order { degree[row + 1] += degree[row] }
        var next = Array(degree[0..<order])
        var neighbours = [Int](repeating: 0, count: degree[order])
        row = 0
        while row < order {
            var position = starts[row] + 1
            let positionEnd = starts[row + 1]
            while position < positionEnd {
                let column = columns[position]
                neighbours[next[row]] = column
                next[row] += 1
                neighbours[next[column]] = row
                next[column] += 1
                position += 1
            }
            row += 1
        }
        return (degree, neighbours)
    }
}
