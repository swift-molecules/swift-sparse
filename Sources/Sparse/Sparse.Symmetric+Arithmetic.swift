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
        for row in 0..<order { result[row] = 0 }
        for row in 0..<order {
            let first = starts[row]
            let x = vector[row]
            var sum = values[first] * x
            for position in (first + 1)..<starts[row + 1] {
                let column = columns[position]
                let value = values[position]
                sum += value * vector[column]
                result[column] += value * x
            }
            result[row] += sum
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
        for row in 0..<order {
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
        }
        starts[order] = columns.count
        return Self(pattern: Sparse.Pattern(order: order, starts: starts, columns: columns), values: sums)
    }
}
