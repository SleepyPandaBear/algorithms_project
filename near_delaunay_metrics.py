import draw
import utils

# NOTE(miha): First metric on quads.
# TODO(miha): decide how we want to return angle, should we return sum of non
# LD angles, list of all non LD angles?
def opposing_angles(points, triangles):
    common_edges = utils.get_common_edges(triangles)

    for e,t in common_edges.items():
        # NOTE(miha): We need to have 2 triangles and they need to be
        # different.
        if len(t) == 2 and t[0] != t[1]:

            # NOTE(miha): New common edge is created of opposing points in a
            # quad.
            new_common_edge = list(set(t[0]) ^ set(t[1]))

            a0 = e[0]
            b0 = new_common_edge[0]
            c0 = e[1]
            angle0 = utils.get_angle(points, a0, b0, c0)

            a1 = e[0]
            b1 = new_common_edge[1]
            c1 = e[1]

            angle1 = utils.get_angle(points, c1, b1, a1)

            #if angle0 + angle1 < 180:
            #    print("not LD")

            print(angle0, angle1, angle0 + angle1)

# NOTE(miha): Second metric on quads.
# TODO(miha): What should we do with ratio? sum it up? put it in a list?
def dual_edge_ratio(points, triangles):
    common_edges = utils.get_common_edges(triangles)

    for e,t in common_edges.items():
        # NOTE(miha): Check if we have a quadliteral
        if len(t) == 2 and t[0] != t[1]:

            # NOTE(miha): Point not in the triangle t[0]
            inscribed = set(t[1]) - set(t[0])
            inscribed = inscribed.pop()

            center, radius = utils.circumcenter(points, t[0])
            in_circle = utils.distance(center, points[inscribed]) <= radius

            # NOTE(miha): If we have non LD edge.
            if in_circle:
                other_inscribed = set(t[0]) - set(t[1])
                other_inscribed = other_inscribed.pop()

                other_center, other_radius = utils.circumcenter(points, t[1])

                uv = utils.distance(points[e[0]], points[e[1]])
                pq = utils.distance(center, other_center)

                print("1:", utils.p_distance(points[e[0]], points[e[1]], 1))
                print("2:", pq/uv)
                print("3:", utils.p_distance(points[e[0]], points[e[1]], 3))
                print("4:", utils.p_distance(points[e[0]], points[e[1]], 4))
                print("5:", utils.p_distance(points[e[0]], points[e[1]], 5))
                print("6:", utils.p_distance(points[e[0]], points[e[1]], 6))
                print("7:", utils.p_distance(points[e[0]], points[e[1]], 7))
                print("70:", utils.p_distance(points[e[0]], points[e[1]], 70))
                print("700:", utils.p_distance(points[e[0]], points[e[1]], 700))
                print("max:", utils.max_distance(points[e[0]], points[e[1]]))
                print("+++++++++++++++++++++++++++++")

# NOTE(miha): Third metric on quads.
def dual_area_overlap(points, triangles):
    common_edges = utils.get_common_edges(triangles)

    for e,t in common_edges.items():
        if len(t) == 2 and t[0] != t[1]:

            # NOTE(miha): Point not in the triangle t[0]
            inscribed = set(t[1]) - set(t[0])
            inscribed = inscribed.pop()

            center, radius = utils.circumcenter(points, t[0])
            in_circle = utils.distance(center, points[inscribed]) <= radius
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

                # NOTE(miha): There is no intersection beteen lines, so we
                # can't calculate its area. (prob we should calculate line0 and
                # line1 lerp at higher number, but happens rarelly)
                if p0 is None or p1 is None:
                    break

                pts = [center, p0, other_center, p1]
                area = utils.calculate_quad_area(*pts)
                print("area:", area)

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

    dual_edge_ratio(points, dt_triangles)
    print("---------")
    dual_edge_ratio(points, non_dt_triangles)

    #dual_area_overlap(points, dt_triangles)
    #print("---------")
    #dual_area_overlap(points, non_dt_triangles)

    #draw.triangles(points, dt_triangles)
    #draw.show()
    #draw.triangles(points, non_dt_triangles)
    #draw.show()

main()
