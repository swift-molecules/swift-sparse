extension Sparse {
    /// A renumbering of a matrix's rows and columns, chosen so that its factor fills in little.
    public struct Ordering: Hashable, Sendable {
        /// `permutation[k]` is the original index placed at position `k`.
        public let permutation: [Int]
        /// `position[i]` is the position of original index `i`: the inverse of `permutation`.
        public let position: [Int]

        @usableFromInline
        internal init(permutation: [Int]) {
            var position = [Int](repeating: 0, count: permutation.count)
            for (index, original) in permutation.enumerated() { position[original] = index }
            self.permutation = permutation
            self.position = position
        }
    }
}

extension Sparse.Ordering {
    /// The identity: every index keeps its place.
    @inlinable
    public static func natural(order: Int) -> Self {
        Self(permutation: Array(0..<order))
    }

    /// A nested dissection of the pattern's graph.
    ///
    /// Each connected piece larger than `leaf` vertices is cut by one level of a breadth-first search from a
    /// pseudo-peripheral vertex, the smallest level near the middle, trimmed of vertices that touch only one side.
    /// The two sides are numbered first, each dissected the same way, and the cut last, so eliminating a side never
    /// fills the other. Pieces of at most `leaf` vertices are numbered in reverse Cuthill–McKee order.
    public static func dissection(of pattern: Sparse.Pattern, leaf: Int = 48) -> Self {
        var dissection = Dissection(pattern: pattern, leaf: Swift.max(leaf, 1))
        return dissection.ordering()
    }
}

extension Sparse.Ordering {
    @usableFromInline
    internal struct Dissection {
        let order: Int
        let starts: [Int]
        let neighbours: [Int]
        let leaf: Int
        /// The piece a vertex currently belongs to.
        var piece: [Int]
        var level: [Int]
        var permutation: [Int]
        var queue: [Int]
        var stamp = 0

        init(pattern: Sparse.Pattern, leaf: Int) {
            let adjacency = pattern.adjacency()
            self.order = pattern.order
            self.starts = adjacency.starts
            self.neighbours = adjacency.neighbours
            self.leaf = leaf
            self.piece = [Int](repeating: 0, count: pattern.order)
            self.level = [Int](repeating: -1, count: pattern.order)
            self.permutation = [Int](repeating: -1, count: pattern.order)
            self.queue = [Int](repeating: 0, count: pattern.order)
        }

        mutating func ordering() -> Sparse.Ordering {
            var work: [(vertices: [Int], first: Int)] = [(Array(0..<order), 0)]
            while let (vertices, first) = work.popLast() {
                stamp += 1
                for vertex in vertices { piece[vertex] = stamp }
                // WHY: a breadth-first search from the first vertex tells whether the piece is connected; if not, it splits into its components.
                let reached = search(from: vertices[0])
                if reached < vertices.count {
                    var start = first
                    let component = Array(queue[0..<reached])
                    work.append((component, start))
                    start += reached
                    let label = piece[vertices[0]]
                    stamp += 1
                    for vertex in component { piece[vertex] = stamp }
                    for vertex in vertices where piece[vertex] == label {
                        let count = search(from: vertex)
                        let next = Array(queue[0..<count])
                        stamp += 1
                        for member in next { piece[member] = stamp }
                        work.append((next, start))
                        start += count
                    }
                    continue
                }
                if vertices.count <= leaf {
                    number(vertices.count, from: first)
                    continue
                }
                guard let cut = separate(vertices, from: first) else {
                    number(vertices.count, from: first)
                    continue
                }
                work.append(contentsOf: cut)
            }
            return Sparse.Ordering(permutation: permutation)
        }

        /// Numbers the `count` vertices last found by `search` in reverse breadth-first order (reverse Cuthill–McKee).
        mutating func number(_ count: Int, from first: Int) {
            var index = 0
            while index < count {
                permutation[first + count - 1 - index] = queue[index]
                index += 1
            }
        }

        /// A breadth-first search within the current piece; fills `queue` and `level`, returns the number reached.
        mutating func search(from root: Int) -> Int {
            let label = piece[root]
            stamp += 1
            let seen = stamp
            piece[root] = seen
            queue[0] = root
            level[root] = 0
            var head = 0
            var tail = 1
            while head < tail {
                let vertex = queue[head]
                head += 1
                var position = starts[vertex]
                let positionEnd = starts[vertex + 1]
                while position < positionEnd {
                    let neighbour = neighbours[position]
                    if piece[neighbour] == label {
                        piece[neighbour] = seen
                        level[neighbour] = level[vertex] + 1
                        queue[tail] = neighbour
                        tail += 1
                    }
                    position += 1
                }
            }
            for index in 0..<tail { piece[queue[index]] = label }
            return tail
        }

        /// Splits a connected piece into two sides and the cut between them, or nil when no level separates it.
        mutating func separate(_ vertices: [Int], from first: Int) -> [(vertices: [Int], first: Int)]? {
            var root = vertices[0]
            var count = search(from: root)
            var depth = level[queue[count - 1]]
            // WHY: a pseudo-peripheral root makes long, thin levels: restart from a vertex of least degree in the last level while the depth grows.
            var round = 0
            while round < 8 {
                var candidate = queue[count - 1]
                var least = Int.max
                var index = count - 1
                while index >= 0, level[queue[index]] == depth {
                    let vertex = queue[index]
                    let degree = starts[vertex + 1] - starts[vertex]
                    if degree < least {
                        least = degree
                        candidate = vertex
                    }
                    index -= 1
                }
                let reached = search(from: candidate)
                let reach = level[queue[reached - 1]]
                if reach <= depth {
                    count = search(from: root)
                    break
                }
                root = candidate
                count = reached
                depth = reach
                round += 1
            }
            guard depth >= 2 else { return nil }
            var sizes = [Int](repeating: 0, count: depth + 1)
            for index in 0..<count { sizes[level[queue[index]]] += 1 }
            // WHY: the cut is the smallest level that leaves between 30% and 70% of the piece before it, else the middle level.
            var chosen = -1
            var middle = -1
            var before = sizes[0]
            var candidate = 1
            while candidate < depth {
                let after = before + sizes[candidate]
                if middle < 0 && 2 * after >= count { middle = candidate }
                if 10 * before >= 3 * count && 10 * after <= 7 * count {
                    if chosen < 0 || sizes[candidate] < sizes[chosen] { chosen = candidate }
                }
                before = after
                candidate += 1
            }
            if chosen < 0 { chosen = middle >= 0 ? middle : depth / 2 }
            let label = piece[root]
            var near: [Int] = []
            var far: [Int] = []
            var cut: [Int] = []
            var index = 0
            while index < count {
                let vertex = queue[index]
                let at = level[vertex]
                if at < chosen {
                    near.append(vertex)
                } else if at > chosen {
                    far.append(vertex)
                } else {
                    // WHY: a vertex of the cut with no neighbour on the far side separates nothing, so it joins the near side.
                    var touches = false
                    var position = starts[vertex]
                    let positionEnd = starts[vertex + 1]
                    while position < positionEnd {
                        let neighbour = neighbours[position]
                        if piece[neighbour] == label && level[neighbour] == chosen + 1 {
                            touches = true
                            break
                        }
                        position += 1
                    }
                    if touches { cut.append(vertex) } else { near.append(vertex) }
                }
                index += 1
            }
            for (index, vertex) in cut.enumerated() {
                permutation[first + near.count + far.count + index] = vertex
            }
            return [(near, first), (far, first + near.count)]
        }
    }
}
