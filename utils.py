from collections import defaultdict
import math
import random
import numpy as np
from pyDelaunay2D.delaunay2D import Delaunay2D

# NOTE(miha): Here are gathered all kinds of functions that we need for our
# project. We have a mixture of math, helper,... functions.

###############################################################################
#                                                                             #
#                             MATH related utils                              #
#                                                                             #
###############################################################################

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
        return (A[0] + (t * Bx_Ax), A[1] + (t * By_Ay)) 

    return None

# NOTE(miha): 'b' is a common point.
def get_angle(points, a, b, c):
    # TODO(miha): Here we need to use different metric aka manhathan, p-norm, infinity
    ang = math.degrees(math.atan2(points[c][1]-points[b][1], points[c][0]-points[b][0]) - math.atan2(points[a][1]-points[b][1], points[a][0]-points[b][0]))
    return ang + 360 if ang < 0 else ang

def circumcenter(points, triangle):
    a = points[triangle[0]]
    b = points[triangle[1]]
    c = points[triangle[2]]

    pts = np.array([a, b, c])

    pts2 = np.dot(pts, pts.T)
    A = np.bmat([[2 * pts2, [[1],
                             [1],
                             [1]]],
                  [[[1, 1, 1, 0]]]])

    b = np.hstack((np.sum(pts * pts, axis=1), [1]))
    x = np.linalg.solve(A, b)
    bary_coords = x[:-1]
    center = np.dot(bary_coords, pts)

    radius = np.linalg.norm(pts[0] - center) # euclidean distance
    #radius = np.sum(np.square(pts[0] - center))  # squared distance
    return (center, radius)

def in_circumcircle(points, triangle, d):
    a = points[triangle[0]]
    b = points[triangle[1]]
    c = points[triangle[2]]

    d = points[d.pop()]

    print("in_circum", a, b, c, d)

    #det = [[a[0]-d[0], a[1]-d[1], (a[0]*a[0]-d[0]*d[0])+(a[1]*a[1]-d[1]*d[1])], \
    #       [b[0]-d[0], b[1]-d[1], (b[0]*b[0]-d[0]*d[0])+(b[1]*b[1]-d[1]*d[1])], \
    #       [c[0]-d[0], c[1]-d[1], (c[0]*c[0]-d[0]*d[0])+(c[1]*c[1]-d[1]*d[1])]]

    det = [[1, a[0], a[1], (a[0]*a[0] + a[1]*a[1])], \
           [1, b[0], b[1], (b[0]*b[0] + b[1]*b[1])], \
           [1, c[0], c[1], (c[0]*c[0] + c[1]*c[1])], \
           [1, d[0], d[1], (d[0]*d[0] + d[1]*d[1])]]

    result = np.linalg.det(det)

    return result

# NOTE(miha): Linear interpolation&exterpolation between two points, returns
# tuple (x,y) which is a point on the interpolated line at 'x'.
def lerp(p0, p1, x):
    y = p0[1] + ((p1[1]-p0[1])/(p1[0]-p0[0]))*(x-p0[0])
    return (x,y)

###############################################################################
#                                                                             #
#                             DELAUNAY related utils                          #
#                                                                             #
###############################################################################

def make_delaunay(points):
    dt = Delaunay2D()
    for p in points:
        dt.addPoint(p)
    return dt

def get_triangles(dt):
    return dt.exportTriangles()

# NOTE(miha): In this function we "skew" Delaunay triangulation, so we can
# measure how close it is to the original DT.
def flip_some_edges(points, triangles, num_of_flips):
    triangles = triangles[:]

    # TODO(miha): We can just keep track of already flipped triangles so we
    # don't flip them again :)

    flips = 0
    flipped_triangles = []
    flipped_edges = []

    for i in range(num_of_flips):
        found_common_edge = False
        tris = triangles[:]

        shuffle0 = triangles[:]
        random.shuffle(triangles)
        shuffle1 = triangles[:]
        random.shuffle(triangles)

        for t0 in shuffle0:
            if found_common_edge:
                break

            for t1 in shuffle1:
                common_edge = list(set(t0) & set(t1))

                # NOTE(miha): If we already flipped this common edge, skip
                # flipping.
                if common_edge in flipped_edges:
                    break

                if t0 != t1 and len(common_edge) == 2                          \
                and t0 not in flipped_triangles and t1 not in flipped_triangles:

                    new_common_edge = list(set(t0) ^ set(t1))

                    intersects = line_intersection(points[common_edge[0]], points[common_edge[1]],
                            points[new_common_edge[0]], points[new_common_edge[1]])

                    # NOTE(miha): If old and new common edge don't intersects,
                    # flip is illegal so skip it.
                    if intersects == None:
                        break

                    tris.remove(t0)
                    tris.remove(t1)

                    new_t0 = (new_common_edge[0], new_common_edge[1], common_edge[0])
                    new_t1 = (new_common_edge[0], new_common_edge[1], common_edge[1])

                    tris.append(new_t0)
                    tris.append(new_t1)

                    flipped_triangles.append(t0)
                    flipped_triangles.append(t1)
                    flipped_edges.append(common_edge)

                    found_common_edge = True

                    flips += 1
                    break

        triangles = tris[:]

    return triangles, flips

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

###############################################################################
#                                                                             #
#                             OTHER related utils                             #
#                                                                             #
###############################################################################

def random_points(n, x_range=(-10, 10), y_range=(-10, 10)):
    points = []

    for _ in range(n):
        points.append((round(random.uniform(*x_range), 2), \
                       round(random.uniform(*y_range), 2)))

    return points

