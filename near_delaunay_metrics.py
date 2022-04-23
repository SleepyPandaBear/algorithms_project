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
        if p0 <= p1:
            edges[(p0, p1)].append((p0,p1,p2))
        else:
            edges[(p1, p0)].append((p0,p1,p2))

        # p1-p2 edge
        if p1 <= p2:
            edges[(p1, p2)].append((p0,p1,p2))
        else:
            edges[(p2, p1)].append((p0,p1,p2))

        # p0-p2 edge
        if p0 <= p2:
            edges[(p0, p2)].append((p0,p1,p2))
        else:
            edges[(p2, p0)].append((p0,p1,p2))

    return edges

def extract_triangles(common_edges):
    triangles = []
    for e,t in common_edges.items():
        triangles.extend(t)

    return triangles

# NOTE(miha): In this function we "skew" Delaunay triangulation, so we can
# measure how close it is to the original DT.
def flip_some_edges(points, triangles, probability=0.1):
    common_edges = get_common_edges(triangles)
    new_common_edges = defaultdict(list)
    print("common edges", common_edges)
    bad_edges = []
    good_edges = []
    flips = 0
    tris = set()

    for e,t in common_edges.items():
        # NOTE(miha): Skip edge edges - can't have common edge with another
        # triangle.
        print("common edges", new_common_edges)
        if len(t) < 2:
            tris.add(t[0])
            new_common_edges[e].extend(t)
            continue

        tris.add(t[0])
        tris.add(t[1])
        
        if random.random() < probability:
            old_edge0 = e[0] 
            old_edge1 = e[1]

            # TODO(miha): Refractor, we can use variable e from for loop!
            new_edge = set(t[0]) ^ set(t[1])
            new_edge0 = new_edge.pop()
            new_edge1 = new_edge.pop()

            if new_edge0 > new_edge1:
                tmp = new_edge1
                new_edge1 = new_edge0
                new_edge0 = tmp

            if line_intersection(points[old_edge0], points[old_edge1], points[new_edge0], points[new_edge1]):
                new_triangle0 = (new_edge0, new_edge1, e[0])
                new_triangle1 = (new_edge0, new_edge1, e[1])
                print("new triangle:", new_triangle0, new_triangle1)
                new_common_edges[e].append(new_triangle0)
                new_common_edges[e].append(new_triangle1)

                flips += 1
        else:
            print("old triangle:", t[0], t[1])
            new_common_edges[e].append(t[0])
            new_common_edges[e].append(t[1])
    
    # NOTE(miha): Add new edges and delete old ones
    for ge, be in zip(good_edges, bad_edges):
        tmp = common_edges[be]
        #new_common_edges[be] = []
        #new_common_edges[ge[0]].append(tmp[0])
        #new_common_edges[ge[0]].append(tmp[1])

    new_triangles = extract_triangles(new_common_edges)
    print("len(t):", len(triangles), " len(nt):", len(new_triangles))
    print("new_common_edge:", new_common_edges)
    return new_triangles

def random_points(n, x_range=(-10, 10), y_range=(-10, 10)):
    points = []

    for _ in range(n):
        points.append((round(random.uniform(*x_range), 2), \
                       round(random.uniform(*y_range), 2)))

    return points

def main():
    points = random_points(5)

    dt = Delaunay2D()
    for p in points:
        dt.addPoint(p)

    dt_triangles = dt.exportTriangles()

    print(points)
    non_dt_triangles = flip_some_edges(points, dt_triangles, probability=1)
    print("dt_triangles:", dt_triangles)
    print("non_dt_triangles:", non_dt_triangles)

    draw_triangles(points, dt_triangles)
    draw_triangles(points, non_dt_triangles)

main()
