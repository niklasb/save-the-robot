# Save the robot

Programming project for the course 'Computational Geometry' at KIT
winter 2015/2016

### Usage

    $ make  # to build C++ program
    $ python3 robot.py <infile> <radius>

### Example

    $ python3 robot.py input/input3.in 80
    Input:
      Radius: 80
      Warehouse: 1038 vertices
      Charging stations: 130

    Computing visibility graph...
    Number of edges in visibility graph: 11068
    Wrote visibility graph to output/visiblity.ipe
    Computing SSSP for interesting points...
    Computing shortest path in reachability graph...
    Wrote reachability graph to output/reachability.ipe
    Goal distance: 3539.8530409000014
    Wrote path to output/path.ipe

    Statistics:
      Time visibility graph: 2.441 sec
      Time reachable graph: 0.023 sec
      Time final path: 0.031 sec

The output can be inspected in [Ipe](http://ipe.otfried.org/) by opening the
file `output/path.ipe`. You can find the output of the above example run
rendered as PDF files in the directory `example_output/`.

### Notes on Ipe input files

The warehouse polygon has to be drawn with `black` stroke and `gray` fill, the
holes are polygons with `black` stroke and `white` fill. The start and end point
as well as the charging stations are disks with colors `green`, `red` and `blue`,
respectively.

Please note that we did not implement translation, so you need to place the
disks and polygons *exactly* without moving them around after.
