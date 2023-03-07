import matplotlib.pyplot as plt
from voronoi import *
from box import *
from model import *
from scipy.spatial import Voronoi
from seeds import *
from hexagon_fct import *


def generate_BDF(directory, r1, mesh_size, t_wall, t_skin, t_total, dim_x, dim_y, mns, bc, x, y):    # Get the points that form the perimeter of the box (from a svg file) (should be a separate function)

    perimeter = box_choice("Box_Clement_V1.svg", dim_x, dim_y)

    # Create Voronoi points to have honeycomb pattern
    filtered_pts = hexagon_pts(r1, dim_x, dim_y, x, y)  # read_points_from_file(path_dots, path_img, dim_x, dim_y)

    # Generate the voronoi cells (the hexagonal walls)
    vor = Voronoi(filtered_pts)

    # Carve a cercle in the middle of the panel
    regions, vertices = voronoi_finite_polygons_2d(vor)              # Gets the regions and vertices from the points (2D)
    polygons_core = get_polygons(regions, vertices, vor, perimeter)  # List of polygons of the voronoi
    polygons = boundary_conditions(polygons_core)                    # Redefines what polygones to keep (inside/outside a cercle or square)

    #verts, faces, faces_sqr = get_faces_and_verts(polygons, t_total-t_skin, t_skin)  # Get verts and faces out of these polygons (for blender/plot graph)

    ## Print faces and verts in a file = Blender allows non-manifold to manifold
    #write_faces_and_verts(verts, faces, faces_sqr)

    # Classify the vertex and segments in a dictionary (faster than the list to retrieve datas)
    Vertex_2D_ID, Segment_2D_ID, Polygon_2D_ID = dict_vertex_2D(polygons)

    # Forme les polygones du panneau en 3D
    Vertex_3D_ID, Segment_3D_ID, Polygon_3D_ID = dict_vertex_3D(Vertex_2D_ID, Segment_2D_ID, Polygon_2D_ID, t_skin, t_total)

    # Build and Show mesh (GMSH)
    n_segments = len(Segment_2D_ID)
    n_polygons = len(Polygon_2D_ID)
    filenameG, filenameE = seeding_V3(Vertex_3D_ID, Segment_3D_ID, Polygon_3D_ID, n_polygons, n_segments, mesh_size, t_total)

    # Générer le bdf
    write_bdf_V2(filenameG, filenameE, directory, t_skin, t_wall, t_total, dim_x, dim_y, mns, bc)  # Séparer les différents aspects en différents

    # Visualiser en 2D le maillage
    # plot_verts(verts, faces, dim_x, dim_y)

    return print("BDF FILE CREATED UNDER " + directory)


#generate_BDF("writing_BDF_Test", 20, 10.0, 0.5, 2.0, 25, 400, 600, None, None, 0, 0)



## Quick infos for modelling
#nozzle = 0.48
#perimeter_len = get_overall_perimeter(verts, faces)
#area = area_of_polygon(perimeter)
#print("Perimeter length : " + str(perimeter_len))
#density = nozzle*(perimeter_len-4*170)/area   # 170 = for square thing
#volume = nozzle*perimeter_len*21+170*170*4    # Test clement
#print("Density : " + str(density))
#print("Volume : " + str(volume))





