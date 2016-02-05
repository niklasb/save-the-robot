import copy
import heapq
import sys
import time

from collections import namedtuple, defaultdict
from lxml import etree

from geo import *
from graph import *

RobotInstance = namedtuple('RobotInstance',
        ('warehouse', 'charging_stations', 'start', 'goal', 'radius'))

def read_xml(fname):
    with open(fname) as f:
        return etree.parse(f)

def write_xml(fname, doc):
    with open(fname, 'wb') as f:
        f.write(etree.tostring(doc))

def parse_poly(path):
    res = []
    for l in path.split('\n'):
        l = l.strip().split()
        if l and l != ['h']:
            p = Point2D(x=int(l[0]), y=int(l[1]))
            # remove duplicate points
            if not res or p != res[-1]:
                res.append(p)
    return SimplePolygon(res)

def parse_ipe(doc):
    holes = []
    stations = []
    outline = None
    for path in doc.xpath('//path[@stroke="black"]'):
        if path.get('fill') == 'white':
            holes += [parse_poly(path.text)]
        else:
            outline = parse_poly(path.text)
    assert outline is not None
    for e in doc.xpath('//use[@name="mark/disk(sx)"]'):
        # we dont support this kind of transform
        assert not e.get('matrix')
        pos = Point2D(*map(int, e.get('pos').split()))
        s = e.get('stroke')
        if s == 'green':
            start = pos
        elif s == 'red':
            goal = pos
        else:
            stations += [pos]
    return RobotInstance(
            warehouse=PolygonWithHoles(outline, holes),
            charging_stations=stations, start=start, goal=goal, radius=radius)

def add_path(doc, path, color, fat=False):
    page = doc.xpath('//page')[0]
    el = etree.Element('path')
    el.set('stroke', color)
    if fat:
        el.set('pen', 'fat')
    el.text = '\n%d %d m\n' % (path[0].x, path[0].y)
    for p in path[1:]:
        el.text += '%d %d l\n' % (p.x, p.y)
    page.append(el)

class Timer():
    def __init__(self):
        self.total = 0

    def start(self):
        self.t0 = time.time()

    def stop(self):
        self.total += time.time() - self.t0


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: %s infile radius' % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    radius = int(sys.argv[2])

    doc = read_xml(infile)
    instance = parse_ipe(doc)

    intermediate_points = instance.charging_stations + [instance.start, instance.goal]
    all_points = list(instance.warehouse.points()) + intermediate_points

    print('Input:')
    print('  Radius: %d' % instance.radius)
    print('  Warehouse: %d vertices' % sum(1 for _ in instance.warehouse.points()))
    print('  Charging stations: %d' % len(instance.charging_stations))
    print()

    timer_visibility = Timer()
    timer_reachability = Timer()
    timer_path = Timer()

    print('Computing visibility graph...')
    timer_visibility.start()
    visibility_graph = compute_visibility_graph_fast(
            instance.warehouse, all_points, instance.radius)
    timer_visibility.stop()
    print('Number of edges in visibility graph: %d' % sum(1 for _ in visibility_graph.edges()))

    doc_visibility = copy.deepcopy(doc)
    drawn = defaultdict(int)
    for a, b, w in visibility_graph.edges():
        canonical = (min(a,b), max(a,b))
        drawn[canonical] += 1
        if drawn[canonical] == 2:
            continue
        add_path(doc_visibility, [a, b], 'red')
    assert all(x == 2 for x in drawn.values())

    write_xml('output/visibility.ipe', doc_visibility)
    print('Wrote visibility graph to output/visiblity.ipe')

    print('Computing SSSP for interesting points...')
    timer_reachability.start()
    dijkstras = {}
    for c in intermediate_points:
        dijkstra = Dijkstra(visibility_graph, c, radius)
        dijkstra.run()
        dijkstras[c] = dijkstra
    timer_reachability.stop()

    print('Computing shortest path in reachability graph...')
    timer_path.start()
    graph_reachable = WeightedGraph(intermediate_points)
    doc_reachability = copy.deepcopy(doc)
    for a, b in each_pair(intermediate_points):
        if dijkstras[a].is_reachable(b):
            d = dijkstras[a].dis[b]
            graph_reachable.add_edge(a, b, d)
            graph_reachable.add_edge(b, a, d)
            add_path(doc_reachability, [a, b], 'red')

    dijkstra = Dijkstra(graph_reachable, instance.start)
    dijkstra.run()
    timer_path.stop()

    write_xml('output/reachability.ipe', doc_reachability)
    print('Wrote reachability graph to output/reachability.ipe')

    if dijkstra.is_reachable(instance.goal):
        print('Goal distance: %s' % dijkstra.dis[instance.goal])
        station_path = dijkstra.path_to(instance.goal)
        cur = instance.start
        complete_path = [instance.start]
        for nxt in station_path:
            complete_path += dijkstras[cur].path_to(nxt)
            cur = nxt

        doc_path = copy.deepcopy(doc)
        add_path(doc_path, complete_path, 'red', fat=True)

        write_xml('output/path.ipe', doc_path)
        print('Wrote path to output/path.ipe')
    else:
        print('Goal not reachable :(')

    print()
    print('Statistics:')
    print('  Time visibility graph: %.3f sec' % timer_visibility.total)
    print('  Time reachable graph: %.3f sec' % timer_reachability.total)
    print('  Time final path: %.3f sec' % timer_path.total)
