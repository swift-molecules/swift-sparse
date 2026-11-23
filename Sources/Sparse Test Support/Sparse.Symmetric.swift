public import Sparse

extension Sparse.Symmetric {
    /// The `order × order` tridiagonal Toeplitz matrix with `diagonal` on the diagonal and `subdiagonal` beside it.
    ///
    /// Its eigenvalues are `diagonal + 2·subdiagonal·cos(pπ/(order+1))` for `p` in `1...order`, with eigenvectors
    /// `sin(p i π/(order+1))`.
    public static func tridiagonal(
        order: Cardinal,
        diagonal: Scalar,
        subdiagonal: Scalar
    ) throws(Sparse.Error) -> Self {
        let order = Int(bitPattern: order)
        var entries: [Entry] = []
        var index = 0
        while index < order {
            entries.append(Self.entry(index, index, diagonal))
            if index + 1 < order { entries.append(Self.entry(index, index + 1, subdiagonal)) }
            index += 1
        }
        return try Self(order: Cardinal(UInt(order)), entries: entries)
    }

    /// The five-point Laplacian on a `rows × columns` grid of interior points with zero boundary values and unit
    /// spacing, numbered row by row: 4 on the diagonal, −1 for each neighbour.
    ///
    /// Its eigenvalues are `μ(p, rows) + μ(q, columns)` with `μ(p, n) = 4 sin²(pπ/(2(n+1)))`.
    public static func laplacian(rows: Cardinal, columns: Cardinal) throws(Sparse.Error) -> Self {
        let rows = Int(bitPattern: rows)
        let columns = Int(bitPattern: columns)
        var entries: [Entry] = []
        var index = 0
        while index < rows * columns {
            let row = index / columns
            let column = index % columns
            entries.append(Self.entry(index, index, 4))
            if column + 1 < columns { entries.append(Self.entry(index, index + 1, -1)) }
            if row + 1 < rows { entries.append(Self.entry(index, index + columns, -1)) }
            index += 1
        }
        return try Self(order: Cardinal(UInt(rows * columns)), entries: entries)
    }

    /// The square of `laplacian(rows:columns:)`: the thirteen-point biharmonic of a simply supported plate.
    ///
    /// Its eigenvalues are `(μ(p, rows) + μ(q, columns))²`.
    public static func biharmonic(rows: Cardinal, columns: Cardinal) throws(Sparse.Error) -> Self {
        let laplacian = try Self.laplacian(rows: rows, columns: columns)
        let order = laplacian.order
        var neighbours = [[(Int, Scalar)]](repeating: [], count: order)
        laplacian.entries.forEach { entry in
            let row = Int(bitPattern: entry.row)
            let column = Int(bitPattern: entry.column)
            neighbours[row].append((column, entry.value))
            if row != column { neighbours[column].append((row, entry.value)) }
        }
        var entries: [Entry] = []
        neighbours.indices.forEach { i in
            neighbours[i].forEach { j, a in
                neighbours[j].filter { $0.0 >= i }.forEach { l, b in entries.append(Self.entry(i, l, a * b)) }
            }
        }
        return try Self(order: Cardinal(UInt(order)), entries: entries)
    }

    /// A positive definite matrix with about `perRow` scattered off-diagonal entries in each row, placed and valued by
    /// a fixed formula of `seed`: values in (−1, 1), and a diagonal that dominates each row by one.
    public static func scattered(order: Cardinal, perRow: Cardinal, seed: UInt64 = 1) throws(Sparse.Error) -> Self {
        let order = Int(bitPattern: order)
        let perRow = Int(bitPattern: perRow)
        var state = seed &* 0x9E37_79B9_7F4A_7C15 | 1
        func next() -> UInt64 {
            state ^= state << 13
            state ^= state >> 7
            state ^= state << 17
            return state
        }
        var entries: [Entry] = []
        var sums = [Scalar](repeating: 0, count: order)
        var draw = 0
        while draw < order * perRow {
            let row = draw / perRow
            let column = Int(next() % UInt64(order))
            let value = Scalar(Double(next() >> 11) * 0x1p-53 * 2 - 1)
            if column != row {
                entries.append(Self.entry(row, column, value))
                sums[row] += Swift.abs(value)
                sums[column] += Swift.abs(value)
            }
            draw += 1
        }
        sums.indices.forEach { row in entries.append(Self.entry(row, row, sums[row] + 1)) }
        return try Self(order: Cardinal(UInt(order)), entries: entries)
    }

    /// Every entry, both triangles, as rows of a dense matrix.
    public var dense: [[Scalar]] {
        var result = [[Scalar]](repeating: [Scalar](repeating: 0, count: order), count: order)
        entries.forEach { entry in
            let row = Int(bitPattern: entry.row)
            let column = Int(bitPattern: entry.column)
            result[row][column] = entry.value
            result[column][row] = entry.value
        }
        return result
    }

    /// The entry adding `value` at `row`, `column`, from machine integers.
    static func entry(_ row: Int, _ column: Int, _ value: Scalar) -> Entry {
        Entry(Ordinal(UInt(row)), Ordinal(UInt(column)), value)
    }
}
