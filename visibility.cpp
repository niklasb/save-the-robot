#include <bits/stdc++.h>
using namespace std;

struct Vec {
    long long x, y;

    long long dot(const Vec& other) const {
        return x * other.x + y * other.y;
    }

    long long cross(const Vec& other) const {
        return x * other.y - y * other.x;
    }

    Vec operator-(const Vec& other) const {
        return Vec{x - other.x, y - other.y};
    }

    Vec operator+(const Vec& other) const {
        return Vec{x + other.x, y + other.y};
    }

    long long norm() const {
        return dot(*this);
    }

    double angle() const {
        return atan2(y, x);
    }

    bool operator<(const Vec& other) const {
        return tie(x, y) < tie(other.x, other.y);
    }

    bool operator==(const Vec& other) const {
        return tie(x, y) == tie(other.x, other.y);
    }

    bool operator!=(const Vec& other) const {
        return !(*this == other);
    }
};

// > 0 for ccw, < 0 for cw, 0 for colinear
long long cross3(const Vec& a, const Vec& b, const Vec& c) {
    return (b - a).cross(c - a);
}

struct Segment {
    Vec a, b;
    long long norm() const {
        return (b - a).norm();
    }
    long long orientation(const Vec& p) const {
        return cross3(a, b, p);
    }
    bool intersect(const Segment& other) const {
        return orientation(other.a) * orientation(other.b) <= 0
            && other.orientation(a) * other.orientation(b) <= 0;
    }
    bool intersect_inner(const Segment& other) const {
        return orientation(other.a) * orientation(other.b) < 0
            && other.orientation(a) * other.orientation(b) < 0;
    }
    bool contains(const Vec& p) const {
        long long proj = (p - a).dot(b - a);
        return orientation(p) == 0 && proj >= 0 && proj <= norm();
    }
    bool operator<(const Segment& other) const {
        return tie(a, b) < tie(other.a, other.b);
    }

    bool operator==(const Segment& other) const {
        return tie(a, b) == tie(other.a, other.b);
    }
};

ostream& operator<<(ostream& o, const Vec& vec) {
    return o << vec.x << "/" << vec.y;
}

ostream& operator<<(ostream& o, const Segment& seg) {
    return o << seg.a << " <-> " << seg.b;
}

struct Polygon {
    vector<Vec> vertices;
    template <typename F>
    void for_edges(F f) const {
        for (size_t i = 0; i < vertices.size(); ++i)
            f(Segment { vertices[i], vertices[(i + 1) % vertices.size()] });
    }
    bool on_border(const Vec& p) const {
        bool is_on_border = false;
        for_edges([&](const Segment& edge) {
            is_on_border |= edge.contains(p);
        });
        return is_on_border;
    }
};

// compares segments crossing a ray by distance to the ray's starting point
struct CmpCrossingSegments {
    // start and direction of ray
    Vec p;
    const Vec& v;
    bool operator()(const Segment& s1, const Segment& s2) {
        Vec v1 = s1.b - s1.a, v2 = s2.b - s2.a;
        auto l = (s1.a - p).cross(v1) * v.cross(v2);
        auto r = (s2.a - p).cross(v2) * v.cross(v1);
        if (l == r) {
            // my lord, this is a mess. we need to make sure that even if
            // at the current sweepline the segments have the same distance,
            // we compare them properly such that their order will remain valid
            // for the next sweepline
            if (s1.a == s2.a)
                return s1.orientation(s2.b) < 0;
            if (s1.b == s2.b)
                return s1.orientation(s1.b) < 0;
            if (s1.b == s2.a)
                return 1;
            if (s2.b == s1.a)
                return 0;
            if (s1.orientation(s2.a) <= 0 && s1.orientation(s2.b) <= 0)
                return 1;
            if (s1.orientation(s2.a) >= 0 && s1.orientation(s2.b) >= 0)
                return 0;
            if (s2.orientation(s1.a) <= 0 && s2.orientation(s1.b) <= 0)
                return 0;
            if (s2.orientation(s1.a) >= 0 && s2.orientation(s1.b) >= 0)
                return 1;
            assert(0);
            return 0;
        }
        return l < r;
    }
};

const int inf = 1e6;

vector<vector<int>> compute_visibility_graph(
        const vector<Polygon>& polys,
        vector<Vec> points,
        long long radius)
{
    int n = points.size();
    map<Vec, int> order;
    for (int i = 0; i < n; ++i)
        order[points[i]] = i;

    vector<vector<int>> graph(n);
    map<Vec, int> rank;
    vector<Vec> orig_points = points;
    for (int ci = 0; ci < n; ++ci) {
        auto center = orig_points[ci];

        points = orig_points;
        swap(points[0], points[ci]);
        sort(begin(points) + 1, end(points), [&](const Vec& a, const Vec& b) {
            return (a - center).angle() < (b - center).angle();
        });
        vector<vector<Vec>> sweeplines(1);
        for (int i = 1; i < n; ++i) {
            rank[points[i]] = sweeplines.size() - 1;
            sweeplines.back().push_back(points[i]);
            if (i + 1 < n && cross3(center, points[i], points[i + 1]))
                sweeplines.emplace_back();
        }

        vector<vector<Segment>> enter(sweeplines.size()), leave(sweeplines.size());
        for (auto& poly: polys) {
            poly.for_edges([&](const Segment& edge) {
                Vec a = edge.a, b = edge.b;
                long long ccw = cross3(center, a, b);
                if (ccw == 0)
                    return;
                if (ccw < 0)
                    swap(a, b);
                assert(cross3(center, a, b) > 0);
                Segment e {a, b};
                int i = rank[a], j = rank[b];
                if (i < j) {
                    enter[i].push_back(e);
                    leave[j].push_back(e);
                } else {
                    enter[0].push_back(e);
                    leave[j].push_back(e);
                    enter[i].push_back(e);
                }
            });
        }

        vector<bool> switch_outsideness(sweeplines.size());
        bool outside = false;
        bool is_on_border = false;
        for (auto& poly: polys)
            is_on_border |= poly.on_border(center);
        if (is_on_border) {
            outside = true;
            Segment ray { center, center + Vec {-inf, -1}};
            for (auto& poly: polys) {
                poly.for_edges([&](const Segment& edge) {
                    if (ray.intersect_inner(edge))
                        outside = !outside;
                    if (edge.contains(center)) {
                        for (const Vec& x: {edge.a, edge.b}) {
                            if (x != center)
                                switch_outsideness[rank[x]] = 1;
                        }
                    }
                });
            }
        }

        Vec cur_sweepline_dir;
        CmpCrossingSegments cmp { center, cur_sweepline_dir };
        set<Segment, CmpCrossingSegments> crossing_edges(cmp);
        for (size_t i = 0; i < sweeplines.size(); ++i) {
            cur_sweepline_dir = sweeplines[i][0] - center;

            for (auto& edge: enter[i])
                crossing_edges.insert(edge);

#ifdef DEBUG
            for (auto it = crossing_edges.begin(); it != crossing_edges.end(); ++it) {
                if (next(it) != crossing_edges.end())  {
                    assert(!cmp(*next(it), *it));
                }
            }
#endif

            bool switched_outsideness = false;
            if (switch_outsideness[i] && outside) {
                outside = false;
                switched_outsideness = true;
            }

            for (auto& p : sweeplines[i]) {
                bool crossed = false;
                if (crossing_edges.size()) {
                    auto edge = *begin(crossing_edges);
                    crossed = cross3(edge.a, p, edge.b) >  0;
                }
                if (!crossed) {
                    long long dis = (p - center).norm();
                    if (dis <= radius * radius && !outside) {
                        graph[order[center]].push_back(order[p]);
                    }
                }
            }

            for (auto& edge: leave[i])
                crossing_edges.erase(edge);
            if (switch_outsideness[i] && !outside && !switched_outsideness)
                outside = true;
        }
    }
    return graph;
}

Vec read_point() {
    Vec res;
    cin >> res.x >> res.y;
    return res;
}

vector<Vec> read_points() {
    int k; cin >> k;
    vector<Vec> res;
    for (int i = 0; i < k; ++i)
        res.push_back(read_point());
    return res;
}

Polygon read_polygon() {
    return Polygon { read_points() };
}

int main() {
    int num_polys = 0;
    cin >> num_polys;
    vector<Polygon> polys;
    for (int i = 0; i < num_polys; ++i)
        polys.push_back(read_polygon());
    // must be a superset of polys
    auto points = read_points();

    long long radius;
    cin >> radius;
    auto graph = compute_visibility_graph(polys, points, radius);
    for (size_t i = 0; i < graph.size(); ++i) {
        for (int j: graph[i])
            cout << i << " " << j << "\n";
    }
}
