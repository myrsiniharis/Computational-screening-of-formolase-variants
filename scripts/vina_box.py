from pymol import cmd, cgo

import sys

center = [float(x) for x in sys.argv[1:4]]
size = [float(x) for x in sys.argv[4:7]]

x0 = center[0] - size[0]/2
x1 = center[0] + size[0]/2
y0 = center[1] - size[1]/2
y1 = center[1] + size[1]/2
z0 = center[2] - size[2]/2
z1 = center[2] + size[2]/2

box = [
    cgo.BEGIN, cgo.LINES,
    cgo.COLOR, 1.0, 0.0, 0.0,

    # bottom
    cgo.VERTEX, x0, y0, z0, cgo.VERTEX, x1, y0, z0,
    cgo.VERTEX, x1, y0, z0, cgo.VERTEX, x1, y1, z0,
    cgo.VERTEX, x1, y1, z0, cgo.VERTEX, x0, y1, z0,
    cgo.VERTEX, x0, y1, z0, cgo.VERTEX, x0, y0, z0,

    # vertical edges
    cgo.VERTEX, x0, y0, z0, cgo.VERTEX, x0, y0, z1,
    cgo.VERTEX, x1, y0, z0, cgo.VERTEX, x1, y0, z1,
    cgo.VERTEX, x1, y1, z0, cgo.VERTEX, x1, y1, z1,
    cgo.VERTEX, x0, y1, z0, cgo.VERTEX, x0, y1, z1,

    # top
    cgo.VERTEX, x0, y0, z1, cgo.VERTEX, x1, y0, z1,
    cgo.VERTEX, x1, y0, z1, cgo.VERTEX, x1, y1, z1,
    cgo.VERTEX, x1, y1, z1, cgo.VERTEX, x0, y1, z1,
    cgo.VERTEX, x0, y1, z1, cgo.VERTEX, x0, y0, z1,

    cgo.END
]

cmd.load_cgo(box, "vina_box")

# to run in PyMOL: pymol -cq vina_box.py -- -8.702, -10.705, 8.418 20 20 20
