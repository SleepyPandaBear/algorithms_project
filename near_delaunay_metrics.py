import math
import random
from pyDelaunay2D.delaunay2D import Delaunay2D
import matplotlib.pyplot as plt
from collections import defaultdict

# TODO(miha): flip some common edges
# TODO(miha): Draw flipped edges
# TODO(miha): Start implementing metrics

def line_intersection(A, B, C, D): 
    Bx_Ax = B[0] - A[0] 
    By_Ay = B[1] - A[1] 
    Dx_Cx = D[0] - C[0] 
    Dy_Cy = D[1] - C[1] 
    
    determinant = (-Dx_Cx * By_Ay + Bx_Ax * Dy_Cy) 
    if abs(determinant) < 1e-20: 
        # special case for parralel lines
        # same Y coordinate
        if A[1] == B[1] == C[1] == D[1]:
            if C[0] <= A[0] <= D[0] or C[0] <= B[0] <= D[0]:
                return A
        return None 

    s = (-By_Ay * (A[0] - C[0]) + Bx_Ax * (A[1] - C[1])) / determinant 
    t = ( Dx_Cx * (A[1] - C[1]) - Dy_Cy * (A[0] - C[0])) / determinant 

    if s >= 0 and s <= 1 and t >= 0 and t <= 1: 
        # TODO(Miha): Do we really need to return location? Is returning true enough?
        return (A[0] + (t * Bx_Ax), A[1] + (t * By_Ay)) 

    return None

# NOTE(miha): triangles are in form of indices to the points TODO
def draw_triangles(points, triangles):
    for p0,p1,p2 in triangles:
        tmp_t = [points[p0], points[p1], points[p2], points[p0]]
        plt.plot(*zip(*tmp_t), color="blue")
    plt.show()

def draw_points(points):
    plt.scatter(*zip(*points))
    plt.show()

def get_common_edges(triangles):
    edges = defaultdict(list)

    # NOTE(miha): We sort edges so we cannot add two orientations of the same
    # edge.
    for p0,p1,p2 in triangles:
        # p0-p1 edge
        edges[tuple(sorted((p0, p1)))].append(tuple(sorted((p0,p1,p2))))

        # p1-p2 edge
        edges[tuple(sorted((p1, p2)))].append(tuple(sorted((p0,p1,p2))))

        # p0-p2 edge
        edges[tuple(sorted((p0, p2)))].append(tuple(sorted((p0,p1,p2))))

    return edges

def extract_triangles(common_edges):
    triangles = []
    for e,t in common_edges.items():
        triangles.extend(t)

    return triangles

def update_common_edge(common_edges, t, ce, t0, t1):
    if t[0] in common_edges[(ce[0], ce[1])]:
        common_edges[(ce[0], ce[1])].remove(t[0])
        common_edges[(ce[0], ce[1])].append(t0)
    if t[1] in common_edges[(ce[0], ce[1])]:
        common_edges[(ce[0], ce[1])].remove(t[1])
        common_edges[(ce[0], ce[1])].append(t1)

def flip_edge(triangles, common_edges):
    new_common_edges = defaultdict(list)

    e,t = random.choice(list(common_edges.items()))

    while(len(t) < 2):
        e,t = random.choice(list(common_edges.items()))

    new_edge = list(set(t[0]) ^ set(t[1]))

    t0 = (new_edge[0], new_edge[1], e[0])
    t1 = (new_edge[0], new_edge[1], e[1])

    # TODO: add to 'new_common_edges'
    # TODO: update all 4 edge edges of quadliteral

    perimeter = list(set(t[0]) | set(t[1]))
    ce0 = list(sorted((perimeter[0], perimeter[1])))
    ce1 = list(sorted((perimeter[1], perimeter[2])))
    ce2 = list(sorted((perimeter[2], perimeter[3])))
    ce3 = list(sorted((perimeter[3], perimeter[0])))

    print("ce:", common_edges)
    print("perimeter", perimeter)
    print("perimeter_new:", list(set(t0) | set(t1)))
    print("ce0", ce0)
    print("ce[ce0]:", common_edges[(ce0[0], ce0[1])])
    print("t:", t)
    print("new t:", t0, t1)

    update_common_edge(common_edges, t, ce0, t0, t1)
    update_common_edge(common_edges, t, ce1, t0, t1)
    update_common_edge(common_edges, t, ce2, t0, t1)
    update_common_edge(common_edges, t, ce3, t0, t1)

    #common_edges[ce0].remove(BAD TRIANGLE)
    #common_edges[ce0].append(GOOD TRIANGLE)


    return common_edges, new_edge

# NOTE(miha): In this function we "skew" Delaunay triangulation, so we can
# measure how close it is to the original DT.
def flip_some_edges(points, triangles, flips, probability=0.1):
    common_edges = get_common_edges(triangles)
    fliped_edges = []

    # iterate over all common edges (common edges has only lists with 2 triangles)
    # we create a new set for storing trinagles
    # if there is prob to flip an edge, we flip it
    # else no
    # we always append sorted new triangles to new set for storing trinagles

    # better idea:
    # make function for one flip
    # argument accepts number of flips, call above function that many times

    for i in range(flips):
        common_edges, fliped_edge = flip_edge(triangles, common_edges)
        print(fliped_edge)
        print(common_edges)
        print("_--------")
        fliped_edges.append(fliped_edge)
    

    new_triangles = extract_triangles(common_edges)
    print(new_triangles)
    return extract_triangles(common_edges)

def random_points(n, x_range=(-10, 10), y_range=(-10, 10)):
    points = []

    for _ in range(n):
        points.append((round(random.uniform(*x_range), 2), \
                       round(random.uniform(*y_range), 2)))

    return points

def main():
    points = random_points(6)

    dt = Delaunay2D()
    for p in points:
        dt.addPoint(p)

    dt_triangles = dt.exportTriangles()

    non_dt_triangles = flip_some_edges(points, dt_triangles, 1, probability=1)

    draw_triangles(points, dt_triangles)
    draw_triangles(points, non_dt_triangles)

main()
