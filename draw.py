import matplotlib.pyplot as plt

# TODO(miha): option to draw multiple plots in the back
# TODO(miha): option to draw to the same plot?

# NOTE(miha): Triangles are given in the form of indices to the point array.
def triangles(points, triangles):
    for p0,p1,p2 in triangles:
        tmp_t = [points[p0], points[p1], points[p2], points[p0]]
        plt.plot(*zip(*tmp_t), color="blue")

def points(points):
    plt.scatter(*zip(*points))

def point(x, y, c):
    plt.scatter(x, y, color=c)

def show():
    plt.show()
