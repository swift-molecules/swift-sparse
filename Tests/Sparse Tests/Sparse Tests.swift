import Foundation
import Sparse
import Sparse_Test_Support
import Testing

enum `Sparse Tests` {}

extension `Sparse Tests` {
    static func dense(_ matrix: [[Double]], times vector: [Double]) -> [Double] {
        matrix.map { row in zip(row, vector).reduce(0) { $0 + $1.0 * $1.1 } }
    }

    static func distance(_ a: [Double], _ b: [Double]) -> Double {
        zip(a, b).reduce(0) { Swift.max($0, Swift.abs($1.0 - $1.1)) }
    }

    static func vector(_ order: Int) -> [Double] {
        (0..<order).map { Double(($0 * 7919) % 101) / 50 - 1 }
    }
}

extension `Sparse Tests` {
    @Suite
    struct Assembly {
        @Test
        func `Entries at the same place add up, from either triangle`() throws {
            let matrix = try Sparse.Symmetric<Double>(
                order: 3,
                entries: [.init(0, 1, 2), .init(1, 0, 3), .init(2, 2, 1), .init(2, 2, 4)]
            )
            #expect(matrix[0, 1] == 5)
            #expect(matrix[1, 0] == 5)
            #expect(matrix[2, 2] == 5)
            #expect(matrix[0, 2] == 0)
            #expect(matrix.diagonal == [0, 0, 5])
            #expect(matrix.pattern.count == 4)
        }

        @Test
        func `A place whose entries cancel stays stored`() throws {
            let matrix = try Sparse.Symmetric<Double>(order: 2, entries: [.init(0, 1, 1), .init(1, 0, -1)])
            #expect(matrix.pattern.columns(in: 0) == [0, 1])
            #expect(matrix[0, 1] == 0)
        }

        @Test
        func `An entry outside the matrix is refused`() {
            #expect(throws: Sparse.Error.index(row: 0, column: 3)) {
                try Sparse.Symmetric<Double>(order: 3, entries: [.init(0, 3, 1)])
            }
            #expect(throws: Sparse.Error.order) {
                try Sparse.Symmetric<Double>(order: -1, entries: [])
            }
        }

        @Test
        func `The stored entries rebuild the same matrix`() throws {
            let matrix = Sparse.Symmetric<Double>.scattered(order: 50, perRow: 4)
            #expect(try Sparse.Symmetric(order: 50, entries: matrix.entries) == matrix)
        }
    }

    @Suite
    struct Arithmetic {
        @Test
        func `The product matches the dense product`() {
            let matrix = Sparse.Symmetric<Double>.scattered(order: 60, perRow: 5, seed: 3)
            let x = `Sparse Tests`.vector(60)
            let expected = `Sparse Tests`.dense(matrix.dense, times: x)
            #expect(`Sparse Tests`.distance(matrix.applied(to: x), expected) < 1e-12)
        }

        @Test
        func `A sum lives on the union of the patterns and adds entry by entry`() throws {
            let a = Sparse.Symmetric<Double>.scattered(order: 40, perRow: 3, seed: 5)
            let b = Sparse.Symmetric<Double>.laplacian(rows: 5, columns: 8)
            let sum = try a.adding(b, scaledBy: -2.5)
            for row in 0..<40 {
                for column in 0..<40 {
                    #expect(sum[row, column] == a[row, column] - 2.5 * b[row, column])
                }
            }
            #expect(try a.adding(b, scaledBy: 0).pattern == sum.pattern)
            #expect(throws: Sparse.Error.order) { try a.adding(.laplacian(rows: 2, columns: 2)) }
        }

        @Test
        func `The norm is the largest absolute row sum`() {
            #expect(Sparse.Symmetric<Double>.laplacian(rows: 4, columns: 4).norm == 8)
            #expect(Sparse.Symmetric<Double>.tridiagonal(order: 5, diagonal: 2, offDiagonal: -1).norm == 4)
        }
    }

    @Suite
    struct Ordering {
        @Test
        func `A dissection is a permutation, also of a graph in several pieces`() {
            var entries = Sparse.Symmetric<Double>.laplacian(rows: 30, columns: 30).entries
            entries += Sparse.Symmetric<Double>.laplacian(rows: 7, columns: 3).entries.map {
                .init($0.row + 900, $0.column + 900, $0.value)
            }
            let matrix = try! Sparse.Symmetric<Double>(order: 925, entries: entries)
            let ordering = Sparse.Ordering.dissection(of: matrix.pattern, leaf: 8)
            #expect(ordering.permutation.sorted() == Array(0..<925))
            for (position, original) in ordering.permutation.enumerated() {
                #expect(ordering.position[original] == position)
            }
        }

        @Test
        func `On a grid a dissection fills in far less than the natural order`() {
            let pattern = Sparse.Symmetric<Double>.laplacian(rows: 100, columns: 100).pattern
            let natural = Sparse.Elimination(pattern, ordering: .natural(order: pattern.order)).count
            let dissected = Sparse.Elimination(pattern).count
            // The natural order fills the band: about 100 entries a column.
            #expect(natural > 900_000)
            #expect(dissected * 3 < natural)
        }
    }

    @Suite
    struct Factorization {
        @Test(arguments: [1, 2, 7, 48, 1000])
        func `A solve inverts the product`(leaf: Int) throws {
            let matrix = Sparse.Symmetric<Double>.scattered(order: 300, perRow: 4, seed: 11)
            let elimination = Sparse.Elimination(matrix.pattern, ordering: .dissection(of: matrix.pattern, leaf: leaf))
            let factorization = try Sparse.Symmetric.Factorization(matrix, elimination: elimination)
            let x = `Sparse Tests`.vector(300)
            #expect(`Sparse Tests`.distance(factorization.solved(matrix.applied(to: x)), x) < 1e-12)
        }

        @Test
        func `An indefinite matrix factors and solves without pivoting where its pivots are not zero`() throws {
            let matrix = try Sparse.Symmetric<Double>.laplacian(rows: 20, columns: 15)
                .adding(.tridiagonal(order: 300, diagonal: 1, offDiagonal: 0), scaledBy: -3.3)
            let factorization = try matrix.factorization()
            let x = `Sparse Tests`.vector(300)
            #expect(`Sparse Tests`.distance(factorization.solved(matrix.applied(to: x)), x) < 1e-9)
        }

        @Test(arguments: [0.5, 1.7, 3.3, 4.0001, 7.9])
        func `The inertia of A − μI counts the eigenvalues below μ`(shift: Double) throws {
            let rows = 13
            let columns = 9
            func mu(_ p: Int, _ n: Int) -> Double {
                let s = sin(Double(p) * Double.pi / Double(2 * (n + 1)))
                return 4 * s * s
            }
            var below = 0
            for p in 1...rows {
                for q in 1...columns where mu(p, rows) + mu(q, columns) < shift { below += 1 }
            }
            let matrix = try Sparse.Symmetric<Double>.laplacian(rows: rows, columns: columns)
                .adding(.tridiagonal(order: rows * columns, diagonal: 1, offDiagonal: 0), scaledBy: -shift)
            let inertia = try matrix.factorization().inertia
            #expect(inertia.negative == below)
            #expect(inertia.positive == rows * columns - below)
        }

        @Test
        func `A zero pivot is refused and named in the matrix's own numbering`() throws {
            let matrix = try Sparse.Symmetric<Double>(order: 3, entries: [.init(0, 0, 1), .init(2, 2, 1)])
            #expect(throws: Sparse.Error.singular(1)) { try matrix.factorization() }
        }

        @Test
        func `A matrix on another pattern is refused`() throws {
            let elimination = Sparse.Elimination(Sparse.Symmetric<Double>.laplacian(rows: 3, columns: 3).pattern)
            #expect(throws: Sparse.Error.pattern) {
                try Sparse.Symmetric.Factorization(
                    Sparse.Symmetric<Double>.tridiagonal(order: 9, diagonal: 2, offDiagonal: -1),
                    elimination: elimination
                )
            }
        }
    }
}
