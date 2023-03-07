# David Lessard - 4 novembre 2021

import os
import subprocess
import matplotlib.pyplot as plt
from writing_bdf import generate_BDF
from model import op2_reading


def model_write_solve(r1, mesh_size=4, t_wall=0.5, t_skin=1.0, t_total=25.0, dim_x=400, dim_y=600, mns=None, bc=None, x=0, y=0):
    # Create directory name
    dirr = "r1-" + str(r1) + "_msh-" + str(mesh_size) + "_tw-" + str(t_wall) + "_ts-" + str(t_skin) + "_tt-" + str(t_total) + "_x-" + str(dim_x) + "_y-" + str(dim_y) + "_nsm-" + str(mns)
    if bc is not None:
        dirr = dirr + "_bc-inplane"
    if x != 0:
        dirr = dirr + "_x-" + str(x)
    if y != 0:
        dirr = dirr + "_y-" + str(y)

    # Create the path
    parent_dir = 'Y:/Documents/NX Simulations/VIBRATION QUASI-ISO ORTHO/'
    directory = os.path.join(parent_dir, dirr)
    if not os.path.isdir(directory):  # Does not redo a test that was already made!
        generate_BDF(directory, r1, mesh_size, t_wall, t_skin, t_total, dim_x, dim_y, mns, bc, x, y)

        # Solve in "directory"
        myBat = open(r'Y:\Documents\GitHub\Python_Nastran\runBatch.bat', 'w+')
        myBat.write('cd "' + directory + '"\n"C:/Program Files/Siemens/Simcenter3D_2020.2/NXNASTRAN/bin/nastranw.exe" "' + directory +'/model_PyNastran.dat" scr=yes old=no delete=f04,log,xdb\n')
        myBat.close()

        # Run the batch file (in cmd prompt)
        subprocess.call("runBatch.bat")
    else:  # if the file already exist
        print("********************************\nThe directory already exist : " + directory + "\n********************************")

    op2_filename = directory + "/model_pynastran.op2"
    return op2_filename

#filename1 = model_write_solve(r1=10, t_skin=1.0, mns=30)
#op2_filename = 'Y:/Documents/NX Simulations/VIBRATION 2/r1-50_msh-7_tw-0.5_ts-2.0_tt-25.0_x-400_y-600_nsm-None/model_pynastran.op2'


# VIBRATION QUASI-ISO ORTHO (folder)
filename1 = model_write_solve(r1=50, mns=30)

plt.show()






