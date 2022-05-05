from collections import defaultdict
import math
import random
import numpy as np
import sys
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

    # NOTE(miha): classic lin alg style for finding angle
    # cos(a) = dot(x,y)/norm(x)*norm(y)

    p0 = points[a]
    p1 = points[b]
    p2 = points[c]

    v0 = np.array([p0[0] - p1[0], p0[1] - p1[1]])
    v1 = np.array([p1[0] - p2[0], p1[1] - p2[1]])

    ang = math.degrees(math.acos(np.dot(v0, v1)/(np.linalg.norm(v0)*np.linalg.norm(v1))))

    return ang

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

def p_distance(p0, p1, p):
    return (math.pow(p0[0] - p1[0], 2) + math.pow(p0[1] - p1[1], 2))**(1/p)

def max_distance(p0, p1):
    return max(abs(p0[0]), abs(p0[1]), abs(p1[0]), abs(p1[1]))

def distance(p0, p1):
    return math.sqrt(math.pow(p0[0] - p1[0], 2) + math.pow(p0[1] - p1[1], 2))

# NOTE(miha): We divide quad into two triangles and then calculate area of both
# triangles and sum it up.
def calculate_quad_area(a, b, c, d):
    # NOTE(miha): projecting a onto b formula: (np.dot(a, b) / np.dot(b, b)) * b

    # NOTE(miha): projecting (a,b) onto (b,d)
    v0 = np.array([a[0]-b[0], a[1]-b[1]])
    v1 = np.array([d[0]-b[0], d[1]-b[1]])
    p0 = (np.dot(v0, v1)/np.dot(v1, v1))*v1

    # NOTE(miha): projecting (c,d) onto (b,d)
    v0 = np.array([c[0]-d[0], c[1]-d[1]])
    v1 = np.array([b[0]-d[0], b[1]-d[1]])
    p1 = (np.dot(v0, v1)/np.dot(v1, v1))*v1
    
    e0 = np.array([p0[0]-a[0], p0[1]-a[1]])
    e1 = np.array([p1[0]-c[0], p1[1]-c[1]])

    # NOTE(miha): l stands for length

    # TODO(miha): Here we should pass a function for calculating distance :)

    e0l = np.linalg.norm(e0)
    e1l = np.linalg.norm(e1)

    common_edge = np.array([d[0]-b[0], d[1]-b[1]])
    common_edgel = np.linalg.norm(common_edge)

    s0 = e0l * (common_edgel/2)
    s1 = e1l * (common_edgel/2)

    return s0 + s1

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

    return triangles, flips, flipped_edges

# NOTE(miha): This function returns a dictionary with keys consisting of edges
# and values consisting of triangles that are adjacent to that edge.
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

# NOTE(miha): Scale distance between points.
def scale_points(points, scale):
    new_points = []

    for p in points:
        x = p[0] * scale
        y = p[1] * scale
        new_points.append((x,y))

    return new_points

# NOTE(miha): This function flips same edges for the same triangles. It is
# usefull when we need to flip same edges on the scaled up version to check if
# the DT metric really works!
def flip_same_edges(points, triangles, flipped_edges):
    triangles = triangles[:]

    shuffle0 = triangles[:]
    shuffle1 = triangles[:]

    flips = 0
    finished = False
    N = len(flipped_edges)

    for t0 in shuffle0:
        if finished:
            break

        for t1 in shuffle1:
            common_edge = list(set(t0) & set(t1))

            if t0 != t1 and len(common_edge) == 2:
                new_common_edge = list(set(t0) ^ set(t1))

                # TODO(miha): Flip should alway be legal right?
                intersects = line_intersection(points[common_edge[0]], points[common_edge[1]],
                        points[new_common_edge[0]], points[new_common_edge[1]])

                # NOTE(miha): If old and new common edge don't intersects,
                # flip is illegal so skip it.
                if intersects == None:
                    break

                if common_edge in flipped_edges:
                    triangles.remove(t0)
                    triangles.remove(t1)

                    new_t0 = (new_common_edge[0], new_common_edge[1], common_edge[0])
                    new_t1 = (new_common_edge[0], new_common_edge[1], common_edge[1])

                    triangles.append(new_t0)
                    triangles.append(new_t1)

                    flips += 1

                    flipped_edges.remove(common_edge)

                    if flips == N:
                        finished = True
                        break

    return triangles

###############################################################################
#                                                                             #
#                             LINE sweep calculation                          #
#                                                                             #
###############################################################################

def distance(back, curr):
    if len(back[curr]) > 1:
        a = back[curr].pop()
        min_dis = math.sqrt(math.pow(curr[0]-a[0],2)+math.pow(curr[1]-a[1],2))
        min_point = a
        back[curr].add(a)
        for b in back[curr]:
            d = math.sqrt(math.pow(curr[0]-b[0],2)+math.pow(curr[1]-b[1],2)) 
            if d < min_dis:
                min_dis = d
                min_point = b
        return min_point
    return None

def triangulate(S, vertical=True):
    T = set()
    
    edges_to_check = set()
    back_points = defaultdict(set)
    
    n = defaultdict(set)
    
    if vertical:
        s = sorted(S)
    else:
        s = sorted(S, reverse=True, key=lambda k: k[1])

    #s = generify(s)
        
    # new point connects to the others
    for i in range(len(s)):
        edge_buffer = set()
        bad_edges = set()
        
        for p in s[0:i]:
            edge_buffer.add((p, s[i]))
            
            for e in edges_to_check:
                if e is not (p,s[i]):
                    li = line_intersection(*e, p, s[i])
                    if li is not None:
                        bad_edges.add((p,s[i]))
                    
        good_edges = set()
        good_edges = edge_buffer - bad_edges

        for p1,p2 in good_edges:
            back_points[s[i]].add(p1) 
        
        # build edges
        edges_to_check |= good_edges  
        
        # build triangles
        num_of_triangles = 0
        l = len(back_points[s[i]])
        print(num_of_triangles, l)
        while num_of_triangles < l-1:
            mp = distance(back_points, s[i])
            for b in back_points[s[i]]:
                if b is mp:
                    continue

                if (b,mp) in edges_to_check:
                    s1,s2,s3 = sorted((b,mp,s[i]))
                    T.add((s1,s2,s3))
                    num_of_triangles += 1
                elif (mp,b) in edges_to_check:
                    s1,s2,s3 = sorted((b,mp,s[i]))
                    T.add((s1,s2,s3))
                    num_of_triangles += 1                
            back_points[s[i]].discard(mp)
        
    return list(edges_to_check), list(T)

###############################################################################
#                                                                             #
#                             OTHER related utils                             #
#                                                                             #
###############################################################################

def random_points(n, x_range=(-10, 10), y_range=(-10, 10), seed=''):
    points = []

    rng = None
    resulting_seed = None

    if seed:
        rng = random.Random(seed)
        resulting_seed = seed
    else:
        seed = random.randrange(sys.maxsize)
        rng = random.Random(seed)
        resulting_seed = seed

    for _ in range(n):
        points.append((round(rng.uniform(*x_range), 2), \
                       round(rng.uniform(*y_range), 2)))

    return points, seed

def random_grid_points(step=1, div=2, x_range=(-10, 10), y_range=(-10, 10), seed=''):
    points = []

    rng = None
    resulting_seed = None

    rng_distance = step/div

    if seed:
        rng = random.Random(seed)
        resulting_seed = seed
    else:
        seed = random.randrange(sys.maxsize)
        rng = random.Random(seed)
        resulting_seed = seed

    width = (x_range[1] - x_range[0])//step
    height = (y_range[1] - y_range[0])//step

    for h in range(height):
        for w in range(width):
            point = (round(h + rng.uniform(-rng_distance, rng_distance), 2), \
                     round(w + rng.uniform(-rng_distance, rng_distance), 2))

            points.append(point)

    return points, seed
