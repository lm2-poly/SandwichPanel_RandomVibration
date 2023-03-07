import os
import numpy as np
import matplotlib.pyplot as plt
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import read_op2
from scipy.signal import find_peaks
import time
import xlsxwriter

def write_bdf_V2(filenameG, filenameE, directory, t_skin, t_wall, t_total, dim_x, dim_y, mns, bc):
    filename3 = "Grid_Element.bdf"    # Classer en QUAD et TRIA (séparer peaux et core)
    grid = []
    f = open(filename3, 'w')
    if bc is None:
        with open(filenameG) as fp:  # Utilise le format du fichier Grid.bdf pour les GRID
            Lines = fp.readlines()
            for i, line in enumerate(Lines):
                if line[0:4] == "GRID":
                    grid.append([float(line[24:32]), float(line[32:40]), float(line[40:48])])
                    f.write(line)
            fp.close()

    else:
        with open(filenameG) as fp:  # Utilise le format du fichier Grid.bdf pour les GRID
            Lines = fp.readlines()
            for line in Lines:
                if line[0:4] == "GRID":
                    # BOUNDARY CONDITIONS
                    if abs(float(line[24:32])) == 200:  # Block x_edges mouvement in the x direction
                        line = line[:-1] + "        1       " + line[-1:]
                    elif abs(float(line[32:40])) == 300:  # Block y_edges mouvement in the y direction
                        line = line[:-1] + "        2       " + line[-1:]
                    grid.append([float(line[24:32]), float(line[32:40]), float(line[40:48])])
                    f.write(line)
            fp.close()

    with open(filenameE) as fp:  # Utiliser le format du fichier Elements.bdf pour les TRIA6 et CQUAD8
        Lines = fp.readlines()
        for i, line in enumerate(Lines):
            if line[0:6] == "CTRIA6":
                angle = orient_ct6(line, grid)
                line = line[0:16] + "2       " + line[24:] + "        " + angle + "\n"  # Tria6 get pshell 2 properties
                f.write(line)
            elif line[0:6] == "CQUAD8":
                line = line[0:16] + "1       " + line[24:] + Lines[i+1][:-1] + "                                0.000000\n"  # CQUAD8 get pshell 1 properties
                f.write(line)
        fp.close()
        f.write("ENDDATA")
    f.close()

    modele = BDF()  # Generate a BDF object (using PyNastran module)
    modele.read_bdf(filename3, xref=False, punch=True)

    # Ajouter le point central (pour relier avec RBE2)
    mid_node_ID = len(modele.nodes) + 1
    modele.add_grid(mid_node_ID, [0., 0., float(t_total/2)])

    nodes_ID_rbe2 = list()
    for key in modele.nodes:  # Identifie les noeuds qui se trouvent sur le périmètre du cercle
        xyz = modele.nodes[key].xyz
        x = float(xyz[0])
        y = float(xyz[1])
        # Trouve les points de
        d = x**2 + y**2
        if d - (330/2)**2 <= 0.5 and key < mid_node_ID:  # (330/2)**2
            nodes_ID_rbe2.append(key)

    modele.add_rbe2(1, mid_node_ID, "123456", nodes_ID_rbe2)
    #for n, node in enumerate(nodes_ID_rbe2):
    #    modele.add_rbe3(n+1, node, "123456", [1.0], [123], mid_node_ID)

    # Ajouter les cartes liées au random analysis
    modele.add_eigrl(1, 20.0, 2000.0, 10, None, 10, None, "MASS")  # Finds 10 first modes
    modele.add_tabdmp1(20, [20.0, 2000.0], [0.02, 0.02])  # Visqueux - Damping 0.02 = 1/Q = 1/50
    modele.add_freq2(101, 20.0, 2000.0, 10)  # défini les points d'évaluation également espacés (logarithmiquement)
    modele.add_freq4(101, 20.0, 2000.0, 0.1, 4)  # défini les points d'évaluation (détermine si la ligne va être lisse ou non)
    modele.add_spc1(1, "123", [mid_node_ID])  # Bloquer le noeud excité (enforced acceleration/SPCD)
    modele.add_rload1(3, 5, 0.0, 0.0, 5, 0.0, 3)  # 10 = RLOAD2_ID, 11 = SPCD_ID, 12 = DABLED1_ID, 3 = Enforced acceleration using SPCD
    falcon_f = [20.0, 100.0, 300.0, 700.0, 800.0, 925.0, 2000.0]  # Falcon Heavy frequency input (linked to the acceleration input)
    falcon_a = [2910.1, 6507.2, 16991.4, 25954.8, 48059.0, 51677.4, 35206.8]  # Falcon Heavy acceleration input
    modele.add_tabled1(5, falcon_f, falcon_a, "LOG", "LOG")  # Input PSD
    modele.add_spcd(5, [mid_node_ID], ["123"], [1.0])  # Enforced acceleration : node, DOFs, scale
    modele.add_spc1(4, "456", [mid_node_ID])  # Bloquer noeud excité (fréquences propres)
    modele.add_dload(10, 1.0, 1.0, [3])  # Référence au rload1
    modele.add_spcadd(101, [1, 4])  # Union de contrainte (pour les fréquences propres et le sous-cas)
    modele.add_cord2r(4, [0.0, 0.0, 12.5], [0.0, 0.0, 13.5], [1.0, 0.0, 12.5])
    modele.add_randps(102, 1, 1, 1., 0., 1)
    modele.add_tabrnd1(1, falcon_f, [0.0044, 0.0044, 0.01, 0.01, 0.03, 0.03, 0.00644], 'LOG', 'LOG')  # [0.291, 0.6507, 1.699, 2.5954, 4.806, 5.168, 3.5207] #Sinon le RMS est 1000 fois trop grand

    # Material properties (e6 = GPa)
    # CORE (walls) = FACUNDO DATA # DIRECTION 1 IS PERPENDICULAR TO THE FIBER, 2 IS ALONG THE FIBERS
    E11 = 15.603e6  # Dans le sens du matériau, perpendiculaire aux FIBRES)
    E22 = 5.686e6  # Perpendiculaire au matériau (dans le sens des FIBRES)
    G12 = 3.163e6  # IN-PLANE SHEAR MODULUS
    G13 = 1.875e6  # SHEAR OUT-OF-PLAN MODULUS (PERPENDICULAIRE AUX FIBRES - PERP AU PLAN)
    G23 = 1.426e6   # SHEAR OUT-OF-PLAN ALIGNÉ AVEC FIBRES (2)
    Nu12 = 0.38  # Ne change pas entre PEEK et PEEK 20 CF (voir dans overleaf je crois)
    rho1 = 1.39e-6  # Densité (kg/mm^3 = g/cm3*e-6)

    mid1 = 1  # PEEK 20CF ORTHOGONAL
    modele.add_mat8(mid1, e11=E11, e22=E22, nu12=Nu12, g12=G12, rho=rho1, g1z=G13, g2z=G23)

    pid1 = 1
    modele.add_pshell(pid1, mid1=mid1, t=t_wall, mid2=mid1, mid3=mid1)

    # SKINS (depending on the parameters - with/out NSM)
    # DIRECTION 1 IS PERPENDICULAR TO THE FIBER, 2 IS ALONG THE FIBERS
    Es = (138/2)*1e6 # Estimation des propriétés quasi-isotropic du laminé
    Gs = 5.7e6
    nu2 = 0.4
    rho2 = 1.32e-6

    mid2 = 2
    modele.add_mat1(mid2, Es, Gs, nu2, rho2)

    pid2 = 2
    if mns is None:
        modele.add_pshell(pid2, mid1=mid2, t=t_skin, mid2=mid2, mid3=mid2, nsm=mns)
    else:
        modele.add_pshell(pid2, mid1=mid2, t=t_skin, mid2=mid2, mid3=mid2, nsm=mns/2/(dim_x*dim_y - np.pi*(330/2)**2))  # (330/2)**2
    # SKINS (quasi-isotropic) = Propriétés de romain
    # E_11 = 138e6
    # E_22 = 10.2e6
    # G_12 = 5.7e6
    # Nu_21 = 0.3
    # rho_1 = 1.32e-6
    #
    # mid2 = 2
    # modele.add_mat8(mid2, e11=E_11, e22=E_22, nu12=Nu_21, g12=G_12, rho=rho_1)
    #
    # pid2 = 2  # Top Skin
    # if mns is None:
    #     if t_skin == 0.5:
    #         modele.add_pcomp(pid2, mids=[mid2, mid2, mid2, mid2], thicknesses=[0.125, 0.125, 0.125, 0.125],
    #                          thetas=[0., 45., -45., 90.], nsm=mns)
    #     if t_skin == 1.0:
    #         modele.add_pcomp(pid2, mids=[mid2, mid2, mid2, mid2, mid2, mid2, mid2, mid2],
    #                          thicknesses=[0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125],
    #                          thetas=[0., 45., -45., 90., 90., -45., 45., 0.], nsm=mns)
    # else:
    #     if t_skin == 0.5:
    #         modele.add_pcomp(pid2, mids=[mid2, mid2, mid2, mid2], thicknesses=[0.125, 0.125, 0.125, 0.125],
    #                          thetas=[0., 45., -45., 90.],
    #                          nsm=mns / 2 / (dim_x * dim_y - np.pi * (330 / 2) ** 2))
    #     if t_skin == 1.0:
    #         modele.add_pcomp(pid2, mids=[mid2, mid2, mid2, mid2, mid2, mid2, mid2, mid2],
    #                          thicknesses=[0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125],
    #                          thetas=[0., 45., -45., 90., 90., -45., 45., 0.],
    #                          nsm=mns / 2 / (dim_x * dim_y - np.pi * (330 / 2) ** 2))



    modele.sol = 111

    print(modele.get_bdf_stats())
    # for property, value in vars(modele).items():  # pour voir les propriétés de l'objets
    #    print(property, ":", value)

    modele.write_bdf('junk.bdf')  # Ajouter SPC, MAT et Propriété

    # pid_to_eids_map = model.get_element_ids_dict_with_pids()
    # CaseControlDeck.get_op2_data()
    # print(op2.get_op2_stats())

    # Reading data from file1
    with open('HeaderDavid.txt') as fp:
        data = fp.read()

    # Reading data from file2
    with open('junk.bdf') as fp:
        data2 = fp.read()

    data += data2

    os.mkdir(directory)
    with open(directory + '/model_PyNastran.dat', 'w') as fp:  # Ajouter le HEADER
        fp.write(data)

def orient_ct6(line, grid):
    g1 = int(line[24:32])-1
    g2 = int(line[32:40])-1
    x1 = grid[g1][0]
    y1 = grid[g1][1]
    x2 = grid[g2][0]
    y2 = grid[g2][1]
    angle = str(((np.arctan2(y2-y1, x2-x1)/np.pi*180 + 90) % 360) - 180) + "0000000"  # NE MARCHE PAS +90 ou -90!!!!! Essayer un modulo -90/90
    angle = angle[0:8]  # Keeps only 8 caracters to put in the NASTRAN card
    return angle

def op2_reading(filepath):
    time.sleep(5)
    # Waits for the .op2 file being solved
    nas_crit = os.path.isfile(filepath[:-3] + "f04")  # Looks it the temporary .f04 file is still there (the file is deleted because of the command in the .bat file, by default .f04 files are note erased)
    while nas_crit:
        time.sleep(5)  # Waits 5 more seconds
        nas_crit = os.path.isfile(filepath[:-3] + "f04")  # Check status again

    # Open .op2 file and extract useful datas
    op2_file = read_op2(filepath, debug=False)  # Open the .op2 file, and create an OP2 object that contains all results information
    freqs = np.array(op2_file.get_result("modal_contribution.cquad8_stress")[(1, 5, 1, 0, 0, '', '')].freqs)  # List of all frequency evaluated
    max_von_core = np.zeros(len(freqs))  # init
    max_von_skin = np.zeros(len(freqs))  # init
    max_von_c_e = np.zeros(len(freqs))  # init
    max_von_s_e = np.zeros(len(freqs))  # init
    max_von = np.zeros(len(freqs))  # init
    max_von_e = np.zeros(len(freqs))  # init
    for i, f in enumerate(freqs):
        # Core max stress [Element-Nodal Stress]
        von_core = list()
        von_core.append(np.max(op2_file.get_result("modal_contribution.cquad8_stress")[(1, 5, 1, 0, 0, '', '')].data[i, 3::10, 7]) / 1000)  # 10 points par éléments
        von_core.append(np.max(op2_file.get_result("modal_contribution.cquad8_stress")[(1, 5, 1, 0, 0, '', '')].data[i, 5::10, 7]) / 1000)  # 10 points par éléments
        von_core.append(np.max(op2_file.get_result("modal_contribution.cquad8_stress")[(1, 5, 1, 0, 0, '', '')].data[i, 7::10, 7]) / 1000)  # 10 points par éléments
        von_core.append(np.max(op2_file.get_result("modal_contribution.cquad8_stress")[(1, 5, 1, 0, 0, '', '')].data[i, 9::10, 7]) / 1000)  # 10 points par éléments
        max_von_core[i] = np.max(von_core)
        # Skin max stress [Element-Nodal Stress]
        von_skin = list()
        von_skin.append(np.max(op2_file.get_result("modal_contribution.ctria6_stress")[(1, 5, 1, 0, 0, '', '')].data[i, 3::8, 7])/1000)  # 8 points par éléments
        von_skin.append(np.max(op2_file.get_result("modal_contribution.ctria6_stress")[(1, 5, 1, 0, 0, '', '')].data[i, 5::8, 7]) / 1000)  # 8 points par éléments
        von_skin.append(np.max(op2_file.get_result("modal_contribution.ctria6_stress")[(1, 5, 1, 0, 0, '', '')].data[i, 7::8, 7]) / 1000)  # 8 points par éléments
        max_von_skin[i] = np.max(von_skin)
        # Global Maximum [Element-Nodal Stress]
        max_von[i] = max(max_von_core[i], max_von_skin[i])
        # Core max stress [Elemental Stress]
        max_von_c_e[i] = np.max(op2_file.get_result("modal_contribution.cquad8_stress")[(1, 5, 1, 0, 0, '', '')].data[i, 1::10, 7])/1000   # 10 points par éléments
        # Skin max stress [Elemental Stress]
        max_von_s_e[i] = np.max(op2_file.get_result("modal_contribution.ctria6_stress")[(1, 5, 1, 0, 0, '', '')].data[i, 1::8, 7])/1000   # 8 points par éléments
        # Global Maximum [Elemental Stress]
        max_von_e[i] = max(max_von_c_e[i], max_von_s_e[i])

    rms_von_center_s = list()
    rms_von_center_c = list()
    rms_von_corner_s = list()
    rms_von_corner_c = list()
    # Evaluate RMS Stress with Elemental (center) then Element-Nodal (corner)
    for i in range(int(op2_file.get_result("rms.ctria6_stress")[2].ntotal/8)):
        s11 = op2_file.get_result("rms.ctria6_stress")[2].data[0, 8 * i + 1, 0] / 1000
        s22 = op2_file.get_result("rms.ctria6_stress")[2].data[0, 8 * i + 1, 1] / 1000
        s12 = op2_file.get_result("rms.ctria6_stress")[2].data[0, 8 * i + 1, 2] / 1000
        rms_von_center_s.append(np.sqrt(s11 ** 2 - s11 * s22 + s22 ** 2 + 3 * s12 ** 2))
        for j in [3, 5, 7]:
            s11 = op2_file.get_result("rms.ctria6_stress")[2].data[0, 8*i+j, 0] / 1000
            s22 = op2_file.get_result("rms.ctria6_stress")[2].data[0, 8*i+j, 1] / 1000
            s12 = op2_file.get_result("rms.ctria6_stress")[2].data[0, 8*i+j, 2] / 1000
            rms_von_corner_s.append(np.sqrt(s11 ** 2 - s11 * s22 + s22 ** 2 + 3 * s12 ** 2))

    for i in range(int(op2_file.get_result("rms.cquad8_stress")[2].ntotal/10)):
        s11 = op2_file.get_result("rms.cquad8_stress")[2].data[0, 10 * i + 1, 0] / 1000
        s22 = op2_file.get_result("rms.cquad8_stress")[2].data[0, 10 * i + 1, 1] / 1000
        s12 = op2_file.get_result("rms.cquad8_stress")[2].data[0, 10 * i + 1, 2] / 1000
        rms_von_center_c.append(np.sqrt(s11 ** 2 - s11 * s22 + s22 ** 2 + 3 * s12 ** 2))
        for j in [3, 5, 7, 9]:
            s11 = op2_file.get_result("rms.cquad8_stress")[2].data[0, 10*i+j, 0] / 1000
            s22 = op2_file.get_result("rms.cquad8_stress")[2].data[0, 10*i+j, 1] / 1000
            s12 = op2_file.get_result("rms.cquad8_stress")[2].data[0, 10*i+j, 2] / 1000
            rms_von_corner_c.append(np.sqrt(s11 ** 2 - s11 * s22 + s22 ** 2 + 3 * s12 ** 2))

    peaks1, _ = find_peaks(max_von, height=0.1)
    peaks_filt = list()
    for p in peaks1:
        if 2 < p < len(freqs)-2:
            if max_von[p] > max_von[p - 2] and max_von[p] > max_von[p + 2]:
                peaks_filt.append(p)
        else:
            peaks_filt.append(p)

    peaks1, _ = find_peaks(max_von_e, height=0.1)
    peaks_filt_e = list()
    for p in peaks1:
        if 2 < p < len(freqs)-2:
            if max_von_e[p] > max_von_e[p - 2] and max_von_e[p] > max_von_e[p + 2]:
                peaks_filt_e.append(p)
        else:
            peaks_filt_e.append(p)

    # Presentation of results (In the consol)
    print("Fréquences des peaks : ")
    print("Von Mises des peaks : ")
    print(np.array([[freqs[peaks_filt]], [max_von[peaks_filt]]]))
    print("Elemental : ")
    print(np.array([[freqs[peaks_filt_e]], [max_von_e[peaks_filt_e]]]))

    print("Maximum RMS Value for an element's corner SKIN :" + str(np.max(rms_von_corner_s)))
    print("Maximum RMS Value for an element's corner CORE :" + str(np.max(rms_von_corner_c)))
    print("Maximum RMS Value for an element SKIN:" + str(np.max(rms_von_center_s)))
    print("Maximum RMS Value for an element CORE:" + str(np.max(rms_von_center_c)))

    # Presentation of results (Exported in the Folder of the simulation)
    plt.figure()
    plt.plot(freqs, max_von)
    plt.title("Contrainte Von Mises Maximum [Element-Nodal]")
    plt.xlabel("Fréquence [Hz]")
    plt.ylabel("Contrainte Von Mises [MPa]")
    plt.plot(freqs[peaks_filt], max_von[peaks_filt], ".k") # Ajoute les points pour montrer les peaks calculés
    plt.savefig(filepath[:-19] + "plot.png")

    plt.figure()
    plt.plot(freqs, max_von_e)
    plt.title("Contrainte Von Mises Maximum [Elemental]")
    plt.xlabel("Fréquence [Hz]")
    plt.ylabel("Contrainte Von Mises [MPa]")
    plt.plot(freqs[peaks_filt_e], max_von_e[peaks_filt_e], ".k") # Ajoute les points pour montrer les peaks calculés
    plt.savefig(filepath[:-19] + "plot_e.png")

    # Save datas in a Excel spreadsheet (for Matlab to reuse)
    with xlsxwriter.Workbook(filepath[:-19] + 'PSD_VonMises.xlsx') as workbook:
        worksheet = workbook.add_worksheet()
        bold = workbook.add_format({'bold': True})

        for row_num, data in enumerate([["Freqs", "Max Von (Skin)", "Max Von (Core)", "Peaks Freqs", "Peaks Skin", "Peaks Core", "RMS Skin", "RMS Core"], freqs,  max_von_skin, max_von_core, freqs[peaks_filt], max_von_skin[peaks_filt], max_von_core[peaks_filt], [np.max(rms_von_corner_s)], [np.max(rms_von_corner_c)], [" "], ["Elemental"], max_von_s_e, max_von_c_e, freqs[peaks_filt_e], max_von_s_e[peaks_filt_e], max_von_c_e[peaks_filt_e], [np.max(rms_von_center_s)], [np.max(rms_von_center_c)]]):
            if 3 < row_num < 7 or 12 < row_num < 16:
                worksheet.write_row(row_num, 0, data, bold)
            else:
                worksheet.write_row(row_num, 0, data)

    return freqs, max_von

# results = read_op2("Y:/Documents/NX Simulations/Nastran_Mesh/r1-3.01_msh-8_tw-0.5_ts-2.0_tt-25.0_x-400_y-600_nsm-None/model_pynastran.op2")