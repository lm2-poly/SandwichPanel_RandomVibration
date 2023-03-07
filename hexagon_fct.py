import numpy as np
import matplotlib.pyplot as plt

def hexagon_pts(r1, dim_x, dim_y, a, b):
    # Regular honeycomb panel, r1 = longueur des côtés
    filtered_pts = list()
    x = list()
    y = list()
    r2 = r1*np.sqrt(3)/2
    n_pts_x = np.linspace(0, int(np.ceil(dim_x / (2 * r2))), int(np.ceil(dim_x / (2 * r2)) + 1))
    n_pts_y = np.linspace(0, int(np.ceil(dim_y / (1.5 * r1))), int(np.ceil(dim_y / (1.5 * r1)) + 1))
    for ny in n_pts_y:
        for nx in n_pts_x:
            filtered_pts.append([-dim_x / 2 + 2 * nx * r2 - a, -dim_y / 2 + 3 * ny * r1 - b])  # rangée 1
            #x.append(-dim_x / 2 + 2 * nx * r2 - a)  # for debug
            #y.append(-dim_y / 2 + 2 * ny * r1 - b)  # for debug
            filtered_pts.append([-dim_x / 2 + (2 * nx + 1) * r2 - a, -dim_y / 2 + (3 * ny + 1.5) * r1 - b])  # rangée 2
            #x.append(-dim_x / 2 + (2 * nx + 1) * r2 - a)  # for debug
            #y.append(-dim_y / 2 + (2 * ny + 1) * r1 - b)  # for debug
    #plt.scatter(x,y)  # for debug
    #plt.show()  # for debug
    return filtered_pts

#hexagon_pts(5, 400, 600, 0, 0)