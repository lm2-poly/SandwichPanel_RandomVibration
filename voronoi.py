# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import cmath
from svgpathtools import svg2paths, wsvg, svg2paths2
from scipy.spatial import Voronoi
from scipy.spatial import distance
from shapely.geometry import MultiPoint, Point, Polygon, LineString
from shapely.validation import explain_validity
from collections import defaultdict


# function to remove Voronoi points outside of the box, otherwise script might crash when calculating Voro
def remove_extra_pts(pts, box):
    filtered_V2 = list()
    for p in pts:
        point = Point(p)
        geom = Polygon(box)
        if point.within(geom):
            filtered_V2.append(p)
    return filtered_V2


# Create voronoi cells - polygon's based
def voronoi_finite_polygons_2d(vor, radius=None):

    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.
    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.
    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.
    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()*2

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


# Get polygons based on regions and vertices. 3rd dimension with box choice.
def get_polygons(regions, vertices, vor, box):

    # Box shape
    box = Polygon(box)

    polygons = []
    for i, region in enumerate(regions):
        polygon = vertices[region]
        poly = Polygon(polygon)
        poly = poly.intersection(box)

        # workaround when class was automatically changed
        if poly.geom_type == 'MultiPolygon':
            for new_polygon in poly:

                polygon = [p for p in new_polygon.exterior.coords]
                polygons.append(polygon)

        else:
            polygon = [p for p in poly.exterior.coords]
            polygons.append(polygon)
    return polygons


def boundary_conditions(poly):
    polygons = list()
    cercle = Point(0, 0).buffer(330/2)  # 330/2 = rayon du cercle¾
    # carre = Polygon([(-70, -70), (70, -70), (70, 70), (-70, 70), (-70, -70)])

    for i, p in enumerate(poly):
        p = Polygon(p)
        # a = p.intersection(cercle) # conserve polygones à l'intérieur de l'anneau
        b = p.difference(cercle)     # conserve polygones à l'extérieur de l'anneau
        # c = b.intersection(carre)
        # d = b.difference(carre)

        #if not a.is_empty:
        #    polygon = [p for p in a.exterior.coords] # for cercle
        #    polygons.append(polygon)
        if not b.is_empty:
            polygon = [p for p in b.exterior.coords]
            polygons.append(polygon)
        # if not c.is_empty:                           # for square
        #    polygon = [p for p in c.exterior.coords]
        #    polygons.append(polygon)
        #if not d.is_empty:
        #    polygon = [p for p in d.exterior.coords]
        #    polygons.append(polygon)


    return polygons


# Function is pretty explicit. Allows to create STL based on verts and faces. Future improvement: dict based (not list)
def get_faces_and_verts(polygons, z_b, z_t):

    coords = list()
    verts = list()
    faces = list()
    faces_sqr = list()


    for polygon in polygons:
        for j, coord in enumerate(polygon):
            X = coord[0]
            Y = coord[1]
            pt = [X, Y]



            #Z_b , Z_t = height_choice(choice, X, Y)
            XYZ_0 = [X, Y, z_b]
            XYZ_1 = [X, Y, z_t]
            pt = [X, Y]


            if pt not in coords:
                coords.append(pt)
                verts.append(XYZ_0)
                verts.append(XYZ_1)

            if j != 0:
                link_a = verts.index(last_coord_0)
                link_b = verts.index(last_coord_1)
                link_c = verts.index(XYZ_1)
                link_d = verts.index(XYZ_0)
                if link_a == link_b and link_c == link_d:
                    print('fuck')
                face_1a = [link_a, link_b, link_c]
                face_1b = [link_d, link_c, link_b]
                #face_1c = [link_b, link_a, link_c]
                #face_1d = [link_b, link_c, link_a]
                #face_1e = [link_c, link_a, link_b]
                #face_1f = [link_c, link_b, link_a]
                face_2 = [link_a, link_c, link_d]
                if face_1a == [2, 3, 1]:
                    a = 5
                if (face_1a not in faces) and (face_1b not in faces): # and (face_1c not in faces) and (face_1d not in faces) and (face_1e not in faces) and (face_1f not in faces):
                    faces.append(face_1a)
                    faces.append(face_2)
                if [link_a, link_b, link_c, link_d] not in faces_sqr and [link_d, link_c, link_b, link_a] not in faces_sqr:
                    faces_sqr.append([link_a, link_b, link_c, link_d])
            last_coord_0 = XYZ_0
            last_coord_1 = XYZ_1

    return verts, faces, faces_sqr


# Alternative - Dict!
def dict_vertex_2D(polygons):

    Vertex_2D_ID = dict()
    Segment_2D_ID = dict()
    Polygon_2D_ID = dict()


    for polygon in polygons:
        for j, coord in enumerate(polygon):
            X = coord[0]  # was round(coord[0], 3)
            Y = coord[1]
            pt = [X, Y]
            if pt not in Vertex_2D_ID.values():
                pos = len(Vertex_2D_ID) + 1
                Vertex_2D_ID[pos] = pt
    Vertex_2D_val = {str(v): k for k, v in Vertex_2D_ID.items()}

    for polygon in polygons:
        for j, coord in enumerate(polygon):
            if j != 0:
                X1 = coord[0]
                Y1 = coord[1]
                pt0 = [X0, Y0]
                pt1 = [X1, Y1]
                if Vertex_2D_val[str(pt0)] < Vertex_2D_val[str(pt1)]:
                    seg = [Vertex_2D_val[str(pt0)], Vertex_2D_val[str(pt1)]]
                else:
                    seg = [Vertex_2D_val[str(pt1)], Vertex_2D_val[str(pt0)]]

                if seg not in Segment_2D_ID.values():

                    pos = len(Segment_2D_ID) + 1
                    Segment_2D_ID[pos] = seg

                X0 = X1
                Y0 = Y1

            else:
                X0 = coord[0]
                Y0 = coord[1]

    Segment_2D_val = {str(v): k for k, v in Segment_2D_ID.items()}

    for i, polygon in enumerate(polygons):
        seg_of_polygon = list()
        for j, coord in enumerate(polygon):

            if j != 0:

                X1 = coord[0]
                Y1 = coord[1]
                pt0 = [X0, Y0]
                pt1 = [X1, Y1]
                if Vertex_2D_val[str(pt0)] < Vertex_2D_val[str(pt1)]:
                    seg = str([Vertex_2D_val[str(pt0)], Vertex_2D_val[str(pt1)]])
                    seg_opposed = 1

                else:
                    seg = str([Vertex_2D_val[str(pt1)], Vertex_2D_val[str(pt0)]])
                    seg_opposed = -1

                value_seg = Segment_2D_val[seg]
                seg_of_polygon.append(seg_opposed*value_seg)


                X0 = X1
                Y0 = Y1

            else:
                X0 = coord[0]
                Y0 = coord[1]


        Polygon_2D_ID[i+1] = seg_of_polygon

    return Vertex_2D_ID, Segment_2D_ID, Polygon_2D_ID





# Represent the rover
def plot_verts(verts, faces, dim_x, dim_y):

    for face in faces:
        dot_0 = face[0]
        dot_1 = face[2]
        x_0 = verts[dot_0][0]
        x_1 = verts[dot_0][1]
        y_0 = verts[dot_1][0]
        y_1 = verts[dot_1][1]

        plt.plot([x_0, y_0], [x_1, y_1], 'k', linewidth=1.5)

        x = [-dim_x/2, dim_x/2, dim_x/2, -dim_x/2]
        y = [-dim_y/2, -dim_y/2, dim_y/2, dim_y/2]
        plt.plot(x, y, 'k-', linewidth=3)
        plt.plot([x[-1], x[0]], [y[-1], y[0]], 'k-', linewidth=3)

    plt.gca().set_aspect('equal')
    plt.show()


# For material use
def get_overall_perimeter(verts, faces):
    d = 0
    for face in faces:
        dot_0 = face[0]
        dot_1 = face[2]
        x_0 = verts[dot_0][0]
        x_1 = verts[dot_0][1]
        y_0 = verts[dot_1][0]
        y_1 = verts[dot_1][1]
        p1 = [x_0, x_1]
        p2 = [y_0, y_1]

        d = d + dist2(p1, p2) / 2  # divided by two bc 2 triangles makes 1 line

    return d


# always useful
def dist2(p1, p2):
    return ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)**0.5


# For stats
def area_of_polygon(box):
    box = Polygon(box)
    return box.area


# This function doesn't has parametric height yet.
def dict_vertex_3D(Vertex_2D_ID, Segment_2D_ID, Polygon_2D_ID, t_skin, t_total):
    Vertex_3D_ID = dict()
    Segment_3D_ID = dict()
    Polygon_3D_ID = dict()

    # Define both skins, then core
    for keys in Vertex_2D_ID:
        Vertex_3D_ID[keys] = Vertex_2D_ID[keys] + [t_skin/2]
    for keys in Vertex_2D_ID:
        Vertex_3D_ID[keys + len(Vertex_2D_ID)] = Vertex_2D_ID[keys] + [t_total - t_skin/2]

    Vertex_3D_val = {str(v): k for k, v in Vertex_3D_ID.items()}

    # for skins
    # bottom
    for keys in Segment_2D_ID:
        Segment_3D_ID[keys] = Segment_2D_ID[keys]
    # top
    for keys in Segment_2D_ID:
        Segment_3D_ID[keys + len(Segment_2D_ID)] = [Segment_2D_ID[keys][0] + len(Vertex_2D_ID), Segment_2D_ID[keys][1] + len(Vertex_2D_ID)]
    # for core
    for keys in Vertex_2D_ID:
        Segment_3D_ID[len(Segment_3D_ID) + 1] = [keys, keys + len(Vertex_2D_ID)]

    Segment_3D_val = {str(v): k for k, v in Segment_3D_ID.items()}

    # for skins
    # bottom
    for keys in Polygon_2D_ID:
        Polygon_3D_ID[keys] = Polygon_2D_ID[keys]
    # top
    for keys in Polygon_2D_ID:
        temp = list()
        values_to_improve = Polygon_2D_ID[keys]
        for v in values_to_improve:
            if v > 0:
                v = v + len(Segment_2D_ID)
            else:
                v = v - len(Segment_2D_ID)
            temp.append(v)

        Polygon_3D_ID[keys + len(Polygon_2D_ID)] = temp

    # for core

    for keys in Segment_2D_ID:
        if True:        # keys not in Segment_not_in_core:
            pt0 = Segment_2D_ID[keys][0]
            pt1 = Segment_2D_ID[keys][1]
            coord0 = Vertex_2D_ID[pt0]
            coord1 = Vertex_2D_ID[pt1]
            #for seg 1:
            c11 = coord1 + [t_skin/2]
            pt11 = Vertex_3D_val[str(c11)]
            c12 = coord1 + [t_total - t_skin/2]
            pt12 = Vertex_3D_val[str(c12)]

            seg_val_1 = [pt11, pt12]
            seg1 = Segment_3D_val[str(seg_val_1)]

            # for seg 3:
            c31 = coord0 + [t_skin/2]
            pt31 = Vertex_3D_val[str(c31)]
            c32 = coord0 + [t_total - t_skin/2]
            pt32 = Vertex_3D_val[str(c32)]

            seg_val_3 = [pt31, pt32]
            seg3 = -Segment_3D_val[str(seg_val_3)]

            seg0 = keys
            seg2 = -(keys + len(Segment_2D_ID))
            poly = [seg0, seg1, seg2, seg3]
            Polygon_3D_ID[len(Polygon_3D_ID) + 1] = poly

    Polygon_3D_val = {str(v): k for k, v in Polygon_3D_ID.items()}

    return Vertex_3D_ID, Segment_3D_ID, Polygon_3D_ID


def find_segment_not_in_core(Vertex_2D_ID, Segment_2D_ID, cercle, carre):
    segment_not_in_core = list()
    for keys in Segment_2D_ID:
        pt0, pt1 = Segment_2D_ID[keys]
        coord0 = Vertex_2D_ID[pt0]
        coord1 = Vertex_2D_ID[pt1]
        line = LineString([Point(coord0), Point(coord1)])

        x = line.intersection(cercle.boundary)
        y = line.intersection(carre.boundary)
        # TOPOLOGICAL PRB HERE! intersection: prends les segments ayant un point en inters. Et le cercle a des prb d'arrondis
        if not x.is_empty and x.geom_type != 'Point':# and not x.geom_type is 'Point':
            segment_not_in_core.append(keys)
        elif x.geom_type == 'Point' and Point(coord0).distance(cercle.boundary) < 1e-8 and Point(coord1).distance(cercle.boundary) < 1e-8:
            segment_not_in_core.append(keys)
        elif not y.is_empty and not y.geom_type == 'Point':
            segment_not_in_core.append(keys)

    # visual
    #for keys in segment_not_in_core:
    #    pt0, pt1 = Segment_2D_ID[keys]
    #    coord0 = Vertex_2D_ID[pt0]
    #    coord1 = Vertex_2D_ID[pt1]
    #
    #    plt.plot((coord0[0], coord1[0]), (coord0[1], coord1[1]))   # plot the circle?
    #plt.show()


    return segment_not_in_core