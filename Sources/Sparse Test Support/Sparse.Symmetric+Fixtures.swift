public import Sparse

extension Sparse.Symmetric {
    /// The `order × order` tridiagonal Toeplitz matrix with `diagonal` on the diagonal and `offDiagonal` beside it.
    ///
    /// Its eigenvalues are `diagonal + 2·offDiagonal·cos(pπ/(order+1))` for `p` in `1...order`, with eigenvectors
    /// `sin(p i π/(order+1))`.
    public static func tridiagonal(order: Int, diagonal: Scalar, offDiagonal: Scalar) -> Self {
        var entries: [Entry] = []
        for index in 0..<order {
            entries.append(Entry(index, index, diagonal))
            if index + 1 < order { entries.append(Entry(index, index + 1, offDiagonal)) }
        }
        return try! Self(order: order, entries: entries)
    }

    /// The five-point Laplacian on a `rows × columns` grid of interior points with zero boundary values and unit
    /// spacing, numbered row by row: 4 on the diagonal, −1 for each neighbour.
    ///
    /// Its eigenvalues are `μ(p, rows) + μ(q, columns)` with `μ(p, n) = 4 sin²(pπ/(2(n+1)))`.
    public static func laplacian(rows: Int, columns: Int) -> Self {
        var entries: [Entry] = []
        for row in 0..<rows {
            for column in 0..<columns {
                let index = row * columns + column
                entries.append(Entry(index, index, 4))
                if column + 1 < columns { entries.append(Entry(index, index + 1, -1)) }
                if row + 1 < rows { entries.append(Entry(index, index + columns, -1)) }
            }
        }
        return try! Self(order: rows * columns, entries: entries)
    }

    /// The square of `laplacian(rows:columns:)`: the thirteen-point biharmonic of a simply supported plate.
    ///
    /// Its eigenvalues are `(μ(p, rows) + μ(q, columns))²`.
    public static func biharmonic(rows: Int, columns: Int) -> Self {
        let laplacian = Self.laplacian(rows: rows, columns: columns)
        let order = rows * columns
        var neighbours = [[(Int, Scalar)]](repeating: [], count: order)
        for entry in laplacian.entries {
            neighbours[entry.row].append((entry.column, entry.value))
            if entry.row != entry.column { neighbours[entry.column].append((entry.row, entry.value)) }
        }
        var entries: [Entry] = []
        for i in 0..<order {
            for (j, a) in neighbours[i] {
                for (l, b) in neighbours[j] where l >= i {
                    entries.append(Entry(i, l, a * b))
                }
            }
        }
        return try! Self(order: order, entries: entries)
    }

    /// A positive definite matrix with about `perRow` scattered off-diagonal entries in each row, placed and valued by
    /// a fixed formula of `seed`: values in (−1, 1), and a diagonal that dominates each row by one.
    public static func scattered(order: Int, perRow: Int, seed: UInt64 = 1) -> Self {
        var state = seed &* 0x9E37_79B9_7F4A_7C15 | 1
        func next() -> UInt64 {
            state ^= state << 13
            state ^= state >> 7
            state ^= state << 17
            return state
        }
        var entries: [Entry] = []
        var sums = [Scalar](repeating: 0, count: order)
        for row in 0..<order {
            for _ in 0..<perRow {
                let column = Int(next() % UInt64(order))
                if column == row { continue }
                let value = Scalar(Double(next() >> 11) * 0x1p-53 * 2 - 1)
                entries.append(Entry(row, column, value))
                sums[row] += Swift.abs(value)
                sums[column] += Swift.abs(value)
            }
        }
        for row in 0..<order { entries.append(Entry(row, row, sums[row] + 1)) }
        return try! Self(order: order, entries: entries)
    }

    /// Every entry, both triangles, as rows of a dense matrix.
    public var dense: [[Scalar]] {
        var result = [[Scalar]](repeating: [Scalar](repeating: 0, count: order), count: order)
        for entry in entries {
            result[entry.row][entry.column] = entry.value
            result[entry.column][entry.row] = entry.value
        }
        return result
    }
}
