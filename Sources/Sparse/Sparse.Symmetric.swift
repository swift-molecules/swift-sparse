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
    /// The `order × order` matrix that is the sum of `entries`, the way a finite-element assembly adds element
    /// matrices into a global one.
    ///
    /// Each entry adds its value to the entry at `{row, column}`, which is the same entry as `{column, row}`:
    /// give each off-diagonal pair once, from either triangle. Entries at the same place add up, and the place is
    /// stored even when they sum to zero. Fails with `.index` for a row or column outside
    /// `0..<order`.
    @inlinable
    public init(order: Cardinal, entries: some Swift.Sequence<Entry>) throws(Sparse.Error) {
        let order = Int(bitPattern: order)
        var rows: [Int] = []
        var columns: [Int] = []
        var values: [Scalar] = []
        for entry in entries {
            let row = Int(bitPattern: entry.row)
            let column = Int(bitPattern: entry.column)
            guard row < order, column < order else { throw .index(row: row, column: column) }
            rows.append(Swift.min(row, column))
            columns.append(Swift.max(row, column))
            values.append(entry.value)
        }
        var index = 0
        while index < order {
            rows.append(index)
            columns.append(index)
            values.append(0)
            index += 1
        }
        self = Self.compressed(order: order, rows: rows, columns: columns, values: values)
    }

    /// Sorts triplets by row then column with two counting passes, then adds up the duplicates.
    @usableFromInline
    internal static func compressed(order: Int, rows: [Int], columns: [Int], values: [Scalar]) -> Self {
        let total = rows.count
        var byColumn = [Int](repeating: 0, count: order + 1)
        for column in columns { byColumn[column + 1] += 1 }
        var index = 0
        while index < order {
            byColumn[index + 1] += byColumn[index]
            index += 1
        }
        var orderByColumn = [Int](repeating: 0, count: total)
        var entry = 0
        while entry < total {
            orderByColumn[byColumn[columns[entry]]] = entry
            byColumn[columns[entry]] += 1
            entry += 1
        }
        var byRow = [Int](repeating: 0, count: order + 1)
        for row in rows { byRow[row + 1] += 1 }
        index = 0
        while index < order {
            byRow[index + 1] += byRow[index]
            index += 1
        }
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
        var row = 0
        while row < order {
            starts[row] = stored.count
            var last = -1
            var position = byRow[row]
            let positionEnd = byRow[row + 1]
            while position < positionEnd {
                let entry = sorted[position]
                if columns[entry] == last {
                    sums[sums.endIndex - 1] += values[entry]
                } else {
                    last = columns[entry]
                    stored.append(last)
                    sums.append(values[entry])
                }
                position += 1
            }
            row += 1
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
    public subscript(row: Ordinal, column: Ordinal) -> Scalar {
        let low = Int(bitPattern: Swift.min(row, column))
        let high = Int(bitPattern: Swift.max(row, column))
        guard let position = pattern.position(row: low, column: high) else { return 0 }
        return values[position]
    }

    /// The stored entries of the upper triangle, row by row and left to right, the diagonal included.
    @inlinable
    public var entries: [Entry] {
        var result: [Entry] = []
        result.reserveCapacity(values.count)
        var row = 0
        while row < order {
            var position = pattern.starts[row]
            let positionEnd = pattern.starts[row + 1]
            while position < positionEnd {
                result.append(Entry(Ordinal(UInt(row)), Ordinal(UInt(pattern.columns[position])), values[position]))
                position += 1
            }
            row += 1
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
        var row = 0
        while row < order {
            var position = pattern.starts[row]
            let positionEnd = pattern.starts[row + 1]
            while position < positionEnd {
                let size = Swift.abs(values[position])
                sums[row] += size
                let column = pattern.columns[position]
                if column != row { sums[column] += size }
                position += 1
            }
            row += 1
        }
        return sums.max() ?? 0
    }
}

extension Sparse.Symmetric {
    /// The product of this matrix and the column `vector`, which has `order` elements.
    @inlinable
    public func applied(to vector: [Scalar]) -> [Scalar] {
        var result = [Scalar](repeating: 0, count: order)
        apply(to: vector, into: &result)
        return result
    }

    /// Writes the product of this matrix and `vector` into `result`; both have `order` elements.
    @inlinable
    public func apply(to vector: [Scalar], into result: inout [Scalar]) {
        precondition(vector.count == order && result.count == order, "a vector of the matrix's order")
        let starts = pattern.starts
        let columns = pattern.columns
        var row = 0
        while row < order {
            result[row] = 0
            row += 1
        }
        row = 0
        while row < order {
            let first = starts[row]
            let x = vector[row]
            var sum = values[first] * x
            var position = first + 1
            let positionEnd = starts[row + 1]
            while position < positionEnd {
                let column = columns[position]
                let value = values[position]
                sum += value * vector[column]
                result[column] += value * x
                position += 1
            }
            result[row] += sum
            row += 1
        }
    }

    /// Every entry times `factor`, on the same pattern.
    @inlinable
    public func scaled(by factor: Scalar) -> Self {
        Self(pattern: pattern, values: values.map { $0 * factor })
    }

    /// This matrix plus `factor` times `other`, stored on the union of their patterns, the way `K − σM` is formed
    /// for a shift `σ`. The union is kept even where entries cancel, so the result's pattern depends only on the two
    /// patterns. Fails with `.order` when the orders differ.
    @inlinable
    public func adding(_ other: Self, scaledBy factor: Scalar = 1) throws(Sparse.Error) -> Self {
        guard other.order == order else { throw .order }
        if other.pattern == pattern {
            var sum = values
            for position in sum.indices { sum[position] += factor * other.values[position] }
            return Self(pattern: pattern, values: sum)
        }
        var starts = [Int](repeating: 0, count: order + 1)
        var columns: [Int] = []
        var sums: [Scalar] = []
        let a = pattern
        let b = other.pattern
        var row = 0
        while row < order {
            starts[row] = columns.count
            var p = a.starts[row]
            var q = b.starts[row]
            let pEnd = a.starts[row + 1]
            let qEnd = b.starts[row + 1]
            while p < pEnd || q < qEnd {
                let left = p < pEnd ? a.columns[p] : Int.max
                let right = q < qEnd ? b.columns[q] : Int.max
                if left == right {
                    columns.append(left)
                    sums.append(values[p] + factor * other.values[q])
                    p += 1
                    q += 1
                } else if left < right {
                    columns.append(left)
                    sums.append(values[p])
                    p += 1
                } else {
                    columns.append(right)
                    sums.append(factor * other.values[q])
                    q += 1
                }
            }
            row += 1
        }
        starts[order] = columns.count
        return Self(pattern: Sparse.Pattern(order: order, starts: starts, columns: columns), values: sums)
    }
}
