import collections
import heapq

class WeightedGraph:
    def __init__(self, nodes):
        self.nodes = nodes
        self.adj = collections.defaultdict(list)

    def add_edge(self, a, b, w):
        self.adj[a].append((b, w))

    def edges_for(self, a):
        return self.adj[a]

    def edges(self):
        for a, adj in self.adj.items():
            for b, w in adj:
                yield a, b, w


class Dijkstra:
    INFINITY = float('inf')

    def __init__(self, graph, start, radius=INFINITY):
        self.graph = graph
        self.start = start
        self.dis = collections.defaultdict(lambda: self.__class__.INFINITY)
        self.pred = {}
        self.radius = radius

    def run(self, dbg=False):
        s = self.start
        h = [(0, s)]
        self.dis[s] = 0
        self.pred[s] = s
        vis = set()
        while h:
            d, x = heapq.heappop(h)
            if x in vis:
                continue
            vis.add(x)
            if dbg:
                print('Processing node %s' % (x,))
            for y, w in self.graph.edges_for(x):
                dy = d + w
                if dbg:
                    print(dy, dy)
                if dy <= self.radius and dy < self.dis[y] - 1e-9:
                    if dbg:
                        print('Adding %s with dist %s' % (y, dy))
                    self.dis[y] = dy
                    self.pred[y] = x
                    heapq.heappush(h, (dy, y))

    def is_reachable(self, x):
        return self.dis[x] != self.__class__.INFINITY

    def reachable_nodes(self):
        return list(filter(self.is_reachable, self.graph.nodes))

    def path_to(self, x):
        assert self.is_reachable(x)
        res = []
        while self.pred[x] != x:
            res.append(x)
            x = self.pred[x]
        return res[::-1]
