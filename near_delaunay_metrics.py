import math
import random
from pyDelaunay2D.delaunay2D import Delaunay2D
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

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

def have_common_edge(t0, t1):
    return set(t0) & set(t1)

# NOTE(miha): In this function we "skew" Delaunay triangulation, so we can
# measure how close it is to the original DT.
# CARE(miha): Can happen that we flip same edge two times...
def flip_some_edges(points, triangles, flips):
    triangles = triangles[:]

    # TODO(miha): We can just keep track of already flipped triangles so we
    # don't flip them again :)

    for i in range(flips):
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
                if t0 != t1 and len(common_edge) == 2:
                    tris.remove(t0)
                    tris.remove(t1)
                    new_common_edge = list(set(t0) ^ set(t1))

                    new_t0 = (new_common_edge[0], new_common_edge[1], common_edge[0])
                    new_t1 = (new_common_edge[0], new_common_edge[1], common_edge[1])

                    tris.append(new_t0)
                    tris.append(new_t1)

                    found_common_edge = True
                    break

        triangles = tris[:]
    return triangles

def random_points(n, x_range=(-10, 10), y_range=(-10, 10)):
    points = []

    for _ in range(n):
        points.append((round(random.uniform(*x_range), 2), \
                       round(random.uniform(*y_range), 2)))

    return points

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

def opposing_angles(points, triangles):
    common_edges = get_common_edges(triangles)

    for e,t in common_edges.items():
        if len(t) == 2 and t[0] != t[1]:
            common_edge = list(set(t[0]) & set(t[1]))
            b = list(set(t[0]) ^ set(t[1]))

            a0 = common_edge[0]
            b0 = b[0]
            c0 = common_edge[1]
            angle0 = get_angle(points, a0, b0, c0)

            a1 = common_edge[0]
            b1 = b[1]
            c1 = common_edge[1]
            angle1 = get_angle(points, a1, b1, c1)

            print(angle0, angle1, angle0 + angle1)

def distance(p0, p1):
    return math.sqrt(math.pow(p0[0] - p1[0], 2) + math.pow(p0[1] - p1[1], 2))

def dual_edge_ratio(points, triangles):
    common_edges = get_common_edges(triangles)

    for e,t in common_edges.items():
        # NOTE(miha): Check if we have a quadliteral
        if len(t) == 2 and t[0] != t[1]:

            # NOTE(miha): Point not in the triangle t[0]
            inscribed = set(t[1]) - set(t[0])
            inscribed = inscribed.pop()

            center, radius = circumcenter(points, t[0])
            in_circle = distance(center, points[inscribed]) <= radius
            #in_circle = math.sqrt(math.pow(center[0] - points[inscribed][0], 2) + math.pow(center[1] - points[inscribed][1], 2)) <= radius
            print(in_circle)

            # NOTE(miha): If we have non LD edge
            if in_circle:
                other_inscribed = set(t[0]) - set(t[1])
                other_inscribed = other_inscribed.pop()

                other_center, other_radius = circumcenter(points, t[1])

                uv = distance(points[e[0]], points[e[1]])
                pq = distance(center, other_center)

                print("ratio:", pq/uv)

# NOTE(miha): Edge consist of two points
def edge_lerp(edge, x):
    # lerp = x*(1-t) + yt

    p0 = edge[0]
    p1 = edge[1]

    # x1 and y1 = p0
    y = p0[1] + ((p1[1]-p0[1])/(p1[0]-p0[0]))*(x-p0[0])

    return (x,y)

# NOTE(miha): We divide quad into two triangles and then calculate area of both
# triangles and sum it up.
def calculate_quad_area(a, b, c, d):
    # projecting a onto b formula: (np.dot(a, b) / np.dot(b, b)) * b

    plt.scatter(*a, color='r')
    plt.scatter(*b, color='g')
    plt.scatter(*c, color='b')
    plt.scatter(*d, color='k')

    # projecting (a,b) onto (b,d)
    v0 = np.array([a[0]-b[0], a[1]-b[1]])
    v1 = np.array([d[0]-b[0], d[1]-b[1]])
    p0 = (np.dot(v0, v1)/np.dot(v1, v1))*v1
    print("vectors 1:", v0, v1)

    # projecting (c,d) onto (b,d)
    v0 = np.array([c[0]-d[0], c[1]-d[1]])
    v1 = np.array([b[0]-d[0], b[1]-d[1]])
    p1 = (np.dot(v0, v1)/np.dot(v1, v1))*v1
    
    # TODO(miha): Connect p0 with a and p1 with c

    e0 = np.array([p0[0]-a[0], p0[1]-a[1]])
    e1 = np.array([p1[0]-c[0], p1[1]-c[1]])

    e0l = np.linalg.norm(e0)
    e1l = np.linalg.norm(e1)

    common_edge = np.array([d[0]-b[0], d[1]-b[1]])
    common_edgel = np.linalg.norm(common_edge)

    s0 = e0l * (common_edgel/2)
    s1 = e1l * (common_edgel/2)

    print("area:", s0+s1)

    plt.scatter(p0[0]+b[0], p0[1]+b[1], color='y')
    plt.scatter(p1[0]+d[0], p1[1]+d[1], color='y')
    print("p:", p0, p1)
    

def dual_area_overlap(points, triangles):
    common_edges = get_common_edges(triangles)

    for e,t in common_edges.items():
        if len(t) == 2 and t[0] != t[1]:

            # NOTE(miha): Point not in the triangle t[0]
            inscribed = set(t[1]) - set(t[0])
            inscribed = inscribed.pop()

            center, radius = circumcenter(points, t[0])
            in_circle = distance(center, points[inscribed]) <= radius
            #in_circle = math.sqrt(math.pow(center[0] - points[inscribed][0], 2) + math.pow(center[1] - points[inscribed][1], 2)) <= radius

            # NOTE(miha): If we have non LD edge
            if in_circle:
                other_inscribed = set(t[0]) - set(t[1])
                other_inscribed = other_inscribed.pop()

                other_center, other_radius = circumcenter(points, t[1])

                # e - common edge
                # inscribed - point in trianle t[1]
                # other_inscribed - point in trianle t[0]

                # NOTE(miha): Calculate left side, here are edges from pu, qu
                # from the paper
                e0 = (points[inscribed], points[e[0]])
                e1 = (points[other_inscribed], points[e[0]])

                # NOTE(miha): This is a point in the center of a line on the
                # perimeter, we connect this point to the circumcircle center.
                middle_point0 = edge_lerp(e0, (points[inscribed][0] + points[e[0]][0])/2)
                middle_point1 = edge_lerp(e1, (points[other_inscribed][0] + points[e[0]][0])/2)

                # NOTE(miha): Our lines are not super long so we need to
                # increase its length to find intersection
                e0 = (middle_point0, center)
                line0 = (edge_lerp(e0, -100), edge_lerp(e0, 100))
                e1 = (middle_point1, other_center)
                line1 = (edge_lerp(e1, -100), edge_lerp(e1, 100))

                # NOTE(miha): Found another point in the quad we calculate
                # final area
                p0 = line_intersection(line0[0], line0[1], line1[0], line1[1])

                # NOTE(miha): Calculate right side
                e0 = (points[inscribed], points[e[1]])
                e1 = (points[other_inscribed], points[e[1]])

                middle_point0 = edge_lerp(e0, (points[inscribed][0] + points[e[1]][0])/2)
                middle_point1 = edge_lerp(e1, (points[other_inscribed][0] + points[e[1]][0])/2)

                e0 = (middle_point0, center)
                line0 = (edge_lerp(e0, -100), edge_lerp(e0, 100))
                e1 = (middle_point1, other_center)
                line1 = (edge_lerp(e1, -100), edge_lerp(e1, 100))

                p1 = line_intersection(line0[0], line0[1], line1[0], line1[1])

                pts = [center, p0, other_center, p1]
                calculate_quad_area(*pts)


def main():
    points = random_points(6)

    dt = Delaunay2D()
    for p in points:
        dt.addPoint(p)

    dt_triangles = dt.exportTriangles()

    non_dt_triangles = flip_some_edges(points, dt_triangles, 1)
    print(dt_triangles)
    print(non_dt_triangles)

    #opposing_angles(points, dt_triangles)
    #print("-------")
    #opposing_angles(points, non_dt_triangles)
    #dual_edge_ratio(points, dt_triangles)
    #print("---------")
    #dual_edge_ratio(points, non_dt_triangles)
    dual_area_overlap(points, dt_triangles)
    print("---------")
    dual_area_overlap(points, non_dt_triangles)

    #draw_triangles(points, dt_triangles)
    draw_triangles(points, non_dt_triangles)

main()
