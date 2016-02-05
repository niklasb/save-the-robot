import collections
import functools
import math
import random
import subprocess

import graph

class Point2D(collections.namedtuple('Point2D', ('x', 'y'))):
    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def cross(self, other):
        return self.x * other.y - self.y * other.x

    def __sub__(self, other):
        return Point2D(self.x - other.x, self.y - other.y)

    def __add__(self, other):
        return Point2D(self.x + other.x, self.y + other.y)

    def __mul__(self, t):
        return Point2D(t * self.x, t * self.y)

    def rotate(self, angle):
        return Point2D(
                math.cos(angle) * self.x - math.sin(angle) * self.y,
                math.sin(angle) * self.x + math.cos(angle) * self.y)

    def norm(self):
        return self.dot(self)

    def angle(self):
        return math.atan2(self.y, self.x)


def cross3(a, b, c):
    return (b - a).cross(c - a)


class Segment(collections.namedtuple('Segment', ('a', 'b'))):
    def point_orientation(self, p):
        return cross3(self.a, self.b, p)

    def norm(self):
        return (self.b - self.a).norm()

    def intersect(self, other):
        return self.point_orientation(other.a) * self.point_orientation(other.b) <= 0 \
            and other.point_orientation(self.a) * other.point_orientation(self.b) <= 0

    def intersect_inner(self, other):
        return self.point_orientation(other.a) * self.point_orientation(other.b) < 0 \
            and other.point_orientation(self.a) * other.point_orientation(self.b) < 0

    def contains_point(self, p):
        proj = (p - self.a).dot(self.b - self.a)
        return self.point_orientation(p) == 0 and proj >= 0 and proj <= self.norm()

    def mid(self):
        d = self.b - self.a
        assert d.x % 2 == 0 and d.y % 2 == 0
        return self.a + Point2D(d.x // 2, d.y // 2)


class SimplePolygon:
    def __init__(self, points):
        for p in points:
            assert p.x % 2 == 0 and p.y % 2 == 0
        self.points = points

    def edges(self):
        n = len(self.points)
        for i in range(n):
            yield Segment(self.points[i], self.points[(i + 1) % n])

    def border_contains_point(self, p):
        for e in self.edges():
            if e.contains_point(p):
                return True
        return False

    def contains_point(self, p):
        if self.border_contains_point(p):
            return True
        yes = 0
        iterations = 5
        for _ in range(iterations):
            inf = 10**6
            ri = random.randint
            ray = Segment(p, Point2D(inf, ri(-inf,inf)))
            cnt = sum(1 for e in self.edges() if e.intersect(ray))
            yes += cnt % 2
        return yes * 2 >= iterations

    def _intersects_edge(self, s):
        for e in self.edges():
            if e.point_orientation(s.a) * e.point_orientation(s.b) < 0 \
                    and s.point_orientation(e.a) * s.point_orientation(e.b) <= 0:
                return True

    # can return False in case s is inside, but contains a polygon vertex
    def segment_inside(self, s):
        return not self._intersects_edge(s) \
            and self.contains_point(s.mid())
            #and not self.border_contains_point(s.mid())

    def segment_outside(self, s):
        return not self._intersects_edge(s) \
            and (not self.contains_point(s.mid()) or self.border_contains_point(s.mid()))


class PolygonWithHoles:
    def __init__(self, outline, holes):
        self.outline = outline
        self.holes = holes

    def points(self):
        yield from self.outline.points
        for h in self.holes:
            yield from h.points

    def edges(self):
        yield from self.outline.edges()
        for h in self.holes:
            yield from h.edges()

    def segment_inside(self, s):
        return self.outline.segment_inside(s) \
            and all(h.segment_outside(s) for h in self.holes)

    def border_contains_point(self, p):
        return any(poly.border_contains_point(p) for poly in [self.outline] + self.holes)


def each_pair(seq):
    n = len(seq)
    for i in range(n):
        for j in range(i + 1, n):
            yield seq[i], seq[j]

def compute_visibility_graph_naive(polygon, extra_points, radius):
    points = list(set(polygon.points()) | set(extra_points))
    g = graph.WeightedGraph(points)
    for p, q in each_pair(points):
        s = Segment(p, q)
        if p != q and s.norm() <= radius**2 and polygon.segment_inside(s):
            w = s.norm()**0.5
            g.add_edge(p, q, w)
            g.add_edge(q, p, w)
    return g

def compute_visibility_graph(polygon, extra_points, radius):
    poly_points = set(polygon.points())
    points = poly_points | set(extra_points)
    g = graph.WeightedGraph(list(points))

    inf = 10**6
    for center in points:
        points.remove(center)

        # sort and group points by angle around center
        circle = sorted(points, key=lambda q: (q - center).angle())
        rank = {}
        n = len(points)
        sweeplines = [[]]
        for i in range(n):
            q = circle[i]
            rank[q] = len(sweeplines) - 1
            sweeplines[-1].append(q)
            if i + 1 != n and cross3(center, q, circle[i+1]):
                sweeplines.append([])

        # prepare edge enter/leave events
        enter, leave = (collections.defaultdict(list) for _ in range(2))
        for edge in polygon.edges():
            a, b = edge.a, edge.b
            ccw = cross3(center, a, b)
            if ccw == 0:
                # colinear, uninteresting edge
                continue
            if ccw < 0:
                a, b = b, a
            e = Segment(a, b)
            # now we know b is ccw w.r.t a
            i, j = rank[a], rank[b]
            if i < j:
                # regular case
                enter[i].append(e)
                leave[j].append(e)
            else:
                # wrap-around
                enter[0].append(e)
                leave[j].append(e)
                enter[i].append(e)

        # prepare outside/inside events in case where center is on the
        # polygon bounday
        switch_outsideness = set()
        outside = False
        if polygon.border_contains_point(center):
            # we need to check whether our initial sweepline is inside or outside the
            # polygon. we use ray casting starting at center with direction (-inf, -1)
            # because its angle is -pi, which is lower than the angle of the first
            # checked p but still negative (wraparound from pi to -pi occurs between
            # 2nd and 3rd quadrant)
            inf = 10**6
            ray = Segment(center, center + Point2D(-inf, -1))
            outside = True
            for edge in polygon.edges():
                if ray.intersect_inner(edge):
                    outside = not outside
                if edge.contains_point(center):
                    a, b = edge.a, edge.b
                    for x in (a, b):
                        if x != center:
                            switch_outsideness.add(rank[x])

        crossing_edges = []
        for i, ps in enumerate(sweeplines):
            for edge in enter[i]:
                crossing_edges.append(edge)
            switched_outsideness = False
            if i in switch_outsideness and outside:
                outside = False
                switched_outsideness = True

            for p in ps:
                for c in crossing_edges:
                    # we know that for all crossing edges, we have
                    # cross3(center, c.a, p) >= 0 and cross3(center, p, c.b) >= 0,
                    # so we only have to check the distance
                    if cross3(c.a, p, c.b) > 0:
                        break
                else:
                    dis = (p - center).norm()
                    if dis <= radius**2 and not outside:
                        g.add_edge(center, p, dis**0.5)

            for edge in leave[i]:
                crossing_edges.remove(edge)
            if i in switch_outsideness and not outside and not switched_outsideness:
                outside = True

        points.add(center)

    return g

def compute_visibility_graph_fast(polygon, extra_points, radius):
    points = list(set(list(polygon.points()) + extra_points))
    inp = []
    polys = [polygon.outline] + polygon.holes
    inp.append(len(polys))
    for poly in polys:
        inp.append(len(poly.points))
        for p in poly.points:
            inp += [p.x, p.y]
    inp.append(len(points))
    for p in points:
        inp += [p.x, p.y]
    inp.append(radius)

    proc = subprocess.Popen(['./visibility'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    inp = ' '.join(map(str, inp)).encode('utf-8')
    #with open('output/inp', 'wb') as f:
        #f.write(inp)
    out, err = proc.communicate(inp)
    g = graph.WeightedGraph(points)
    res = list(map(int, out.split()))
    for i, j in list(zip(res, res[1:]))[::2]:
        g.add_edge(points[i], points[j], (points[j] - points[i]).norm()**0.5)
    return g
