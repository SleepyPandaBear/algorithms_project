import math
import random
from pyDelaunay2D.delaunay2D import Delaunay2D
from collections import defaultdict
import numpy as np
import draw
import utils

# removed extract_triangles
# removed update_common_edge
# removed flip_edge
# removed have_common_edge

def opposing_angles(points, triangles):
    common_edges = utils.get_common_edges(triangles)

    for e,t in common_edges.items():
        if len(t) == 2 and t[0] != t[1]:
            common_edge = list(set(t[0]) & set(t[1]))
            b = list(set(t[0]) ^ set(t[1]))

            a0 = common_edge[0]
            b0 = b[0]
            c0 = common_edge[1]
            angle0 = utils.get_angle(points, a0, b0, c0)

            a1 = common_edge[0]
            b1 = b[1]
            c1 = common_edge[1]
            angle1 = utils.get_angle(points, a1, b1, c1)

            print(angle0, angle1, angle0 + angle1)

def distance(p0, p1):
    return math.sqrt(math.pow(p0[0] - p1[0], 2) + math.pow(p0[1] - p1[1], 2))

def dual_edge_ratio(points, triangles):
    common_edges = utils.get_common_edges(triangles)

    for e,t in common_edges.items():
        # NOTE(miha): Check if we have a quadliteral
        if len(t) == 2 and t[0] != t[1]:

            # NOTE(miha): Point not in the triangle t[0]
            inscribed = set(t[1]) - set(t[0])
            inscribed = inscribed.pop()

            center, radius = utils.circumcenter(points, t[0])
            in_circle = distance(center, points[inscribed]) <= radius
            #in_circle = math.sqrt(math.pow(center[0] - points[inscribed][0], 2) + math.pow(center[1] - points[inscribed][1], 2)) <= radius
            print(in_circle)

            # NOTE(miha): If we have non LD edge
            if in_circle:
                other_inscribed = set(t[0]) - set(t[1])
                other_inscribed = other_inscribed.pop()

                other_center, other_radius = utils.circumcenter(points, t[1])

                uv = distance(points[e[0]], points[e[1]])
                pq = distance(center, other_center)

                print("ratio:", pq/uv)

# NOTE(miha): We divide quad into two triangles and then calculate area of both
# triangles and sum it up.
def calculate_quad_area(a, b, c, d):
    # projecting a onto b formula: (np.dot(a, b) / np.dot(b, b)) * b

    draw.point(*a, 'r')
    draw.point(*b, 'g')
    draw.point(*c, 'b')
    draw.point(*d, 'k')

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

    draw.point(p0[0]+b[0], p0[1]+b[1], 'y')
    draw.point(p1[0]+d[0], p1[1]+d[1], 'y')
    print("p:", p0, p1)
    

def dual_area_overlap(points, triangles):
    common_edges = utils.get_common_edges(triangles)

    for e,t in common_edges.items():
        if len(t) == 2 and t[0] != t[1]:

            # NOTE(miha): Point not in the triangle t[0]
            inscribed = set(t[1]) - set(t[0])
            inscribed = inscribed.pop()

            center, radius = utils.circumcenter(points, t[0])
            in_circle = distance(center, points[inscribed]) <= radius
            #in_circle = math.sqrt(math.pow(center[0] - points[inscribed][0], 2) + math.pow(center[1] - points[inscribed][1], 2)) <= radius

            # NOTE(miha): If we have non LD edge
            if in_circle:
                other_inscribed = set(t[0]) - set(t[1])
                other_inscribed = other_inscribed.pop()

                other_center, other_radius = utils.circumcenter(points, t[1])

                # e - common edge
                # inscribed - point in trianle t[1]
                # other_inscribed - point in trianle t[0]

                # NOTE(miha): Calculate left side, here are edges from pu, qu
                # from the paper
                e0 = (points[inscribed], points[e[0]])
                e1 = (points[other_inscribed], points[e[0]])

                # NOTE(miha): This is a point in the center of a line on the
                # perimeter, we connect this point to the circumcircle center.
                middle_point0 = utils.lerp(*e0, (points[inscribed][0] + points[e[0]][0])/2)
                middle_point1 = utils.lerp(*e1, (points[other_inscribed][0] + points[e[0]][0])/2)

                # NOTE(miha): Our lines are not super long so we need to
                # increase its length to find intersection
                e0 = (middle_point0, center)
                line0 = (utils.lerp(*e0, -100), utils.lerp(*e0, 100))
                e1 = (middle_point1, other_center)
                line1 = (utils.lerp(*e1, -100), utils.lerp(*e1, 100))

                # NOTE(miha): Found another point in the quad we calculate
                # final area
                p0 = utils.line_intersection(line0[0], line0[1], line1[0], line1[1])

                # NOTE(miha): Calculate right side
                e0 = (points[inscribed], points[e[1]])
                e1 = (points[other_inscribed], points[e[1]])

                middle_point0 = utils.lerp(*e0, (points[inscribed][0] + points[e[1]][0])/2)
                middle_point1 = utils.lerp(*e1, (points[other_inscribed][0] + points[e[1]][0])/2)

                e0 = (middle_point0, center)
                line0 = (utils.lerp(*e0, -100), utils.lerp(*e0, 100))
                e1 = (middle_point1, other_center)
                line1 = (utils.lerp(*e1, -100), utils.lerp(*e1, 100))

                p1 = utils.line_intersection(line0[0], line0[1], line1[0], line1[1])

                pts = [center, p0, other_center, p1]
                calculate_quad_area(*pts)


def main():
    points = utils.random_points(6)

    dt = utils.make_delaunay(points)
    dt_triangles = utils.get_triangles(dt)

    # TODO(miha): Start printing random seeds, so we can get some interesting
    # triangulations later!
    non_dt_triangles, _ = utils.flip_some_edges(points, dt_triangles, 3)
    print(dt_triangles)
    print(non_dt_triangles)

    #opposing_angles(points, dt_triangles)
    #print("-------")
    #opposing_angles(points, non_dt_triangles)
    #dual_edge_ratio(points, dt_triangles)
    #print("---------")
    #dual_edge_ratio(points, non_dt_triangles)
    #dual_area_overlap(points, dt_triangles)
    #print("---------")
    #dual_area_overlap(points, non_dt_triangles)

    draw.triangles(points, dt_triangles)
    draw.show()
    draw.triangles(points, non_dt_triangles)
    draw.show()

main()
