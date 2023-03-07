# coding=utf-8
import scipy
import calfem.geometry as cfg
import calfem.mesh as cfm
import calfem.vis as cfv
from math import cos, sin, pi

from numpy import ceil
import gmsh


def seeding_V3(Vertex_3D_ID, Segment_3D_ID, Polygon_3D_ID, n_polygons, n_segments, mesh_size, t_total):

    gmsh.initialize()
    #gmsh.option.setNumber("General.Terminal", 1)  # Should the information be printed on the terminal (default = 0)

    gmsh.model.add("modele_test")

    # Importer la géométrie dans gmsh.model

    for keys in Vertex_3D_ID:
        gmsh.model.geo.addPoint(Vertex_3D_ID[keys][0], Vertex_3D_ID[keys][1], Vertex_3D_ID[keys][2], 5, keys)

    for keys in Segment_3D_ID:
        gmsh.model.geo.addLine(Segment_3D_ID[keys][0], Segment_3D_ID[keys][1], keys)

    for keys in Polygon_3D_ID:
        gmsh.model.geo.addCurveLoop(Polygon_3D_ID[keys], keys)

    for i, keys in enumerate(Polygon_3D_ID):
        gmsh.model.geo.addPlaneSurface([keys], keys)

    # THIS IS FOR THE STRUCTURED GRID TO APPEAR (VERTICAL ELEMENTS)
    for i, keys in enumerate(Segment_3D_ID):
        if keys > 2 * n_segments:
            gmsh.model.geo.mesh.setTransfiniteCurve(keys, round(t_total/mesh_size)+2)  # Divisions verticales du mesh
    # CQUAD FORMATION, FROM TRIAD
    for i, keys in enumerate(Polygon_3D_ID):
        if keys > 2 * n_polygons:
            gmsh.model.geo.mesh.setTransfiniteSurface(keys)  # Homogeinity of vertical 2D mesh
            gmsh.model.geo.mesh.setRecombine(2, keys)  # Crée des CQUAD en fusionnant des triangles

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.Smoothing", 10)

    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1) # CQUAD8 instead of CQUAD9 (because mesh is 2nd order)

    # ##### OptimizeMesh : Optimize the current mesh with the given algorithm (currently "Gmsh" for default tetrahedral
    # mesh optimizer, "Netgen" for Netgen optimizer, "HighOrder" for direct high-order mesh optimizer,
    # "HighOrderElastic" for high-order elastic smoother, "HighOrderFastCurving" for fast curving algorithm, "Laplace2D"
    # for Laplace smoothing, "Relocate2D" and "Relocate3D" for node relocation).

    # Mesh.QuadqsTopologyOptimizationMethods # 0: default (all),100: pattern-based CAD faces,010: disk quadrangulation
    # remeshing,001: cavity remeshing,xxx: combination of multiple methods (e.g. 111 for all)

    # gmsh.option.setNumber("Mesh.Algorithm", 11)  # Type de mesh des skins : 11=Quad ish

    gmsh.option.setNumber('Mesh.SurfaceFaces', 1)  # Ver las "caras" de los elementos finitos 2D
    gmsh.option.setNumber('Mesh.Points', 1)  # Crée une noeud à chaque place que je crois

    gmsh.option.setNumber('Mesh.MeshSizeMin', mesh_size)
    gmsh.option.setNumber('Mesh.MeshSizeMax', mesh_size)

    gmsh.option.setNumber('Mesh.ElementOrder', 2)  # Set mesh order (Mesh.ElementOrder)

    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.Format", 31)  # Output a BDF

    # Y finalmente guardar la malla
    filenameG = 'Grid.bdf'
    gmsh.write(filenameG)

    gmsh.option.setNumber("Mesh.BdfFieldFormat", 2)  # =2 : Enl<eve les +E12323 à la fin d'une ligne lorsqu'elle continue sur la ligne après, mais fait que c'est des grids*
    filenameE = 'Element.bdf'
    gmsh.write(filenameE)
    # GMSH VISUALISATION (plot of the meshing)
    # gmsh.fltk.run()

    gmsh.finalize()

    return filenameG, filenameE

