# coding=utf-8
import os
import numpy as np
from svgpathtools import svg2paths

# Function that extracts the outline (points and lines) of the geometry from an SVG file to a list
def box_choice(filename, dim_x, dim_y):
    # open svg
    path = os.path.dirname(os.path.abspath(__file__)) + "\Box_SVG\\" + filename
    box = list()
    paths, attributes = svg2paths(path)

    x_final = list()
    y_final = list()

    stop = 50           # debugging
    pts = 20            # 20 pts per arc

    for i, p in enumerate(paths):
        x = list()
        y = list()
        if i < stop:
            # get segments
            for li in p._segments:
                # if line sequence, simply get line
                if type(p._segments[0]).__name__ == "Line":
                    coord = li.start
                    # real represents x, imag y. From svg package.
                    xi = coord.real
                    yi = coord.imag
                    x.append(xi)
                    y.append(yi)
                # complex chapes are defined with CubicBezier. Quick function to represent it
                elif type(p._segments[0]).__name__ == "CubicBezier":
                    X0 = li[0].real
                    X1 = li[1].real
                    X2 = li[2].real
                    X3 = li[3].real
                    Y0 = li[0].imag
                    Y1 = li[1].imag
                    Y2 = li[2].imag
                    Y3 = li[3].imag
                    for j in range(pts + 1):
                        t = j / pts
                        x_bez = (1 - t) ** 3 * X0 + 3 * (1 - t) ** 2 * t * X1 + 3 * (1 - t) * t ** 2 * X2 + (
                                t ** 3) * X3
                        y_bez = (1 - t) ** 3 * Y0 + 3 * (1 - t) ** 2 * t * Y1 + 3 * (1 - t) * t ** 2 * Y2 + (
                                t ** 3) * Y3
                        x.append(x_bez)
                        y.append(y_bez)

        # Write first term
        if i == 0:
            x_final.extend(x)
            y_final.extend(y)
        # Check length of first vs last pt of the segment = allows to reverse if in wrong side
        else:
            distance_premier = ((x[0] - x_last[-1]) ** 2 + (y[0] - y_last[-1]) ** 2) ** 0.5
            distance_dernier = ((x[-1] - x_last[-1]) ** 2 + (y[-1] - y_last[-1]) ** 2) ** 0.5

            # Reverse if wrong side
            if distance_premier > distance_dernier:
                x.reverse()
                y.reverse()

            x_final.extend(x)
            y_final.extend(y)

        x_last = x
        y_last = y

    # get max values for scaling
    x_max = max(x_final)
    x_min = min(x_final)
    y_max = max(y_final)
    y_min = min(y_final)
    delta_x = x_max - x_min
    delta_y = y_max - y_min
    for k in range(len(x_final)):
        # centering, scaling
        x_new = round((x_final[k] - (x_max + x_min) / 2) * dim_x / delta_x)
        y_new = round((y_final[k] - (y_max + y_min) / 2) * dim_y / delta_y)
        coord = [x_new, y_new]

        if coord not in box:
            box.append(coord)

    return box


# old fashion way: define your box point by point
def box_choice_old(choice):
    if choice == 'TestCarreCisaillement':
        box = [(-85, -85), (-85, 85), (85, 85), (85, -85)]

    if choice == 'RoverV2':
        path = 'D:\\Documents\\Uni\\PhD\\000-Git\\GitHub\\OD_PhD_LM2\\venv\\Box\\RoverV2_contourV3_Python.svg'
        box = list()
        paths, attributes = svg2paths(path)

        x_final = list()
        y_final = list()

        stop = 50
        pts = 20

        for i, p in enumerate(paths):
            x = list()
            y = list()
            if i < stop:
                for li in p._segments:
                    if type(p._segments[0]).__name__ == "Line":
                        coord = li.start
                        xi = coord.real
                        yi = coord.imag
                        x.append(xi)
                        y.append(yi)
                    elif type(p._segments[0]).__name__ == "CubicBezier":
                        X0 = li[0].real
                        X1 = li[1].real
                        X2 = li[2].real
                        X3 = li[3].real
                        Y0 = li[0].imag
                        Y1 = li[1].imag
                        Y2 = li[2].imag
                        Y3 = li[3].imag
                        for j in range(pts + 1):
                            t = j / pts
                            x_bez = (1 - t) ** 3 * X0 + 3 * (1 - t) ** 2 * t * X1 + 3 * (1 - t) * t ** 2 * X2 + (
                                        t ** 3) * X3
                            y_bez = (1 - t) ** 3 * Y0 + 3 * (1 - t) ** 2 * t * Y1 + 3 * (1 - t) * t ** 2 * Y2 + (
                                        t ** 3) * Y3
                            x.append(x_bez)
                            y.append(y_bez)

            if i == 0:
                x_final.extend(x)
                y_final.extend(y)
            else:
                distance_premier = ((x[0] - x_last[-1]) ** 2 + (y[0] - y_last[-1]) ** 2) ** 0.5
                distance_dernier = ((x[-1] - x_last[-1]) ** 2 + (y[-1] - y_last[-1]) ** 2) ** 0.5

                if distance_premier > distance_dernier:
                    x.reverse()
                    y.reverse()

                x_final.extend(x)
                y_final.extend(y)

            x_last = x
            y_last = y


        x_max = max(x_final)
        x_min = min(x_final)
        y_max = max(y_final)
        y_min = min(y_final)
        for k in range(len(x_final)-1):
            x_new = round(x_final[k] - (x_max + x_min) / 2)
            y_new = round(y_final[k] - (y_max + y_min) / 2)
            coord = [x_new/4, y_new/4]
            if coord not in box:
                box.append(coord)






    if choice == 'board_Cristina_OD':
        longueur = 1320
        largeur = 220
        theta = 180
        dx = longueur / 2
        dy = largeur / 2
        box = list()

        for i in range(theta+1):
            x = dx + dy*np.sin(np.radians(i))
            y = -dy*np.cos(np.radians(i))
            box.append([x, y])

        for i in range(theta+1):
            x = -dx - dy*np.sin(np.radians(i))
            y = +dy*np.cos(np.radians(i))
            box.append([x, y])

    if choice == 'scutoids':
        box = [(-5, -5), (-5, 5), (5, 5), (5, -5)]

    if choice == 'board_Cristina':
        longueur = 630
        largeur = 220
        theta = 180
        a = 60  # degre 1'
        b = 180 - a  # degre 2'
        dx = longueur / 2
        dy = largeur / 2
        d = 21  # 20.1 #20.05 #19.55  #38 #largeur partie plate en y
        box = list()

        dc = 120  # distance debut cutout
        rayon = 350  # rayon cutout
        x = dx - 180
        y = -dy
        box.append([x, y])

        # cutout1 bas droite
        for i in range(-35, -10, 1):
            x = dx + rayon * np.sin(np.radians(i)) + dc
            y = + rayon * np.cos(np.radians(i)) - 3.6 * dy
            box.append([x, y])
            # x = +dx + dy * np.sin(np.radians(i))
            # y = -dy * np.cos(np.radians(i))
            # box.append([x, y])

        # boutbasdroite
        cx2 = x
        a1 = cx2
        cy2 = y + 20
        a2 = -cy2
        rayon2 = 20
        for i in range(180):
            x = cx2 + rayon2 * np.sin(np.radians(i))
            y = cy2 - rayon2 * np.cos(np.radians(i))
            box.append([x, y])

        # y = - d / 2
        # box.append([x, y])

        c = 42.265  # distance axe sym vertical et partie plate

        e = longueur / 2
        f = 76  # longueur du plat en x
        bb = 21.725  # bloc largeur
        bl = 85  # bloc longeur tot
        bl2 = 31.525  # 31.0525
        bl3 = 32.2  # 32.475 #32.7 #32.975 #32.475
        bl4 = 21.25  # 21.5 #21.75 #22  #21.5
        x1 = e + c
        y1 = - d / 2
        x = x1
        y = y1
        box.append([x, y])  # p1 bas droite

        # y = y1 - bb
        # box.append([x, y])  # pa bas droite
        # xg = x1 - bl2
        # box.append([x, y])  # pg
        # y = y + bb
        # box.append([x, y])  # ph

        # xi = x1 - bl2 - bl3
        # rc= bl2/2
        # px=(x1+xg)/2
        # py= y1
        # for i in range(90,270,1):
        # x = px + rc* np.sin(np.radians(i))
        # y = py + rc * np.cos(np.radians(i))
        # box.append([x, y])

        # xi = x1-bl2-bl3
        # x=xi
        # box.append([x, y])  # pi
        # y = y - bb
        # box.append([x, y])  # pj
        # x = x - bl4
        # box.append([x, y])  # pb bas gauche

        ajoutlong = 0  # 5 #85+ajout long =base triangle
        decalage = 0  # 4.513 #distance base triangle au côté rectangle
        distx1m = 10.804  # 37.039
        disty1m = 30.5  # 30 #47.306
        distx2m = 15.512
        disty2m = disty1m

        xptm = x1 - distx2m
        x = xptm
        yptm = y1 - disty2m
        y = yptm
        box.append([x, y])  # pmtriangle2 bas droite
        x = x1 - bl2
        y = y1
        box.append([x, y])  # p2.3
        y = y1
        x = x1 - bl2 - bl3
        x11 = x
        box.append([x, y])  # p1.1
        xptm = x11 - distx1m
        x = xptm
        yptm = y1 - disty1m
        y = yptm
        box.append([x, y])  # pmtriangle1 bas droite
        x = x1 - bl2 - bl3 - bl4
        x2 = x
        y = y1
        box.append([x, y])  # p2 bas gauche

        # y = d / 2 + decalage
        y3 = d / 2
        y = y3
        x3 = x2
        x = x3
        box.append([x, y])  # p3 haut gauche
        xptm = x3 + distx1m
        x = xptm
        yptm = y3 + disty1m
        y = yptm
        box.append([x, y])  # pmtriangle2
        x = x3 + bl4
        y = y3
        box.append([x, y])  # triangle3
        x = x3 + bl4 + bl3
        x21 = x
        y = y3
        box.append([x, y])  # triangle2.1

        xptm2 = x21 + distx2m
        yptm2 = y3 + disty2m
        x = xptm2
        y = yptm2
        box.append([x, y])  # triangle2m
        # y = y + bb
        # box.append([x, y])  # pc haut gauche
        # x = x +bl4
        # box.append([x, y])  # pu
        # y = y - bb
        # box.append([x, y])  # pv
        # x = x + bl3
        # box.append([x, y])  # pw
        # y = y+ bb
        # box.append([x, y])  # px
        # x = x + bl2
        # box.append([x, y])  # pd haut gauche
        x = x3 + bl4 + bl3 + bl2
        y = y3
        box.append([x, y])  # p4 haut droite

        # bouthautdroite
        for i in range(179):
            x = a1 + rayon2 * np.sin(np.radians(i))
            y = a2 - rayon2 * np.cos(np.radians(i))
            box.append([x, y])

        # cutout2 haut droite
        rayon = 350
        for i in range(-11, -36, -1):
            x = dx + rayon * np.sin(np.radians(i)) + dc
            y = - rayon * np.cos(np.radians(i)) + 3.6 * dy
            box.append([x, y])

        x = dx - 180
        y = -(-dy)
        box.append([x, y])
        # for i in range(b, theta+1):
        # x = dx + dy * np.sin(np.radians(i))
        # y = -dy * np.cos(np.radians(i))
        # box.append([x, y])

        # trou pour la main
        debut = 4
        x = debut
        y = dy
        box.append([x, y])  # point de separation

        lin = 20  # distance de descente a partir du bord
        x = debut
        y = dy - lin
        box.append([x, y])  # point trou debut

        # forme du trou (point 2 à fin-1)
        lt = 120  # long partie plate du trou main
        ht = 30  # hauteur du trou
        x = (lt / 2) - debut
        box.append([x, y])  # point debut courbe

        cx = debut + lt / 2  # coordonée en x du centre
        r = ht / 2  # rayon du cercle
        cy = dy - lin - ht / 2  # coordonée en y du centre
        for i in range(theta, 0, -1):
            x = cx + r * np.sin(np.radians(i))
            y = cy - r * np.cos(np.radians(i))
            box.append([x, y])

        # cote miroir trou
        cx = - lt / 2
        for i in range(theta):
            x = cx - r * np.sin(np.radians(i))
            y = cy - r * np.cos(np.radians(i))
            box.append([x, y])

        x = x + lt / 2
        y = dy - lin
        box.append([x, y])  # point trou fin

        y = dy
        box.append([x, y])  # point merge back

        # côté miroir
        # for i in range(a + 1):
        # x = -dx - dy * np.sin(np.radians(i))
        # y = +dy * np.cos(np.radians(i))
        # box.append([x, y])
        x = -(dx - 180)
        y = -(-dy)
        box.append([x, y])
        rayon = 350
        for i in range(145, 170, 1):
            x = -dx + rayon * np.sin(np.radians(i)) - dc
            y = + rayon * np.cos(np.radians(i)) + 3.6 * dy
            box.append([x, y])

        for i in range(180):
            x = -a1 - rayon2 * np.sin(np.radians(i))
            y = a2 + rayon2 * np.cos(np.radians(i))
            box.append([x, y])

        x = -(e + c)
        y = d / 2
        x1g = x
        y = y3
        box.append([x, y])  # p1 haut gauche
        xptm = x1g + distx2m
        x = xptm
        yptm = y3 + disty2m
        y = yptm
        box.append([x, y])  # pmtriangle2 haut gauche
        x = x1g + bl2
        y = y3
        box.append([x, y])  # p2.3
        y = y3
        x = x1g + bl2 + bl3
        x11 = x
        box.append([x, y])  # p1.1
        xptm = x11 + distx1m
        x = xptm
        yptm = y3 + disty1m
        y = yptm
        box.append([x, y])  # pmtriangle1 haut gauche
        x = x1g + bl2 + bl3 + bl4
        x2g = x
        y = y3
        box.append([x, y])  # p2 haut gauche

        y = y1
        x3 = x2g
        x = x3
        box.append([x, y])  # p3 bas gauche
        xptm = x3 - distx1m
        x = xptm
        yptm = y1 - disty1m
        y = yptm
        box.append([x, y])  # pmtriangle1 bas droite
        x = x3 - bl4
        y = y1
        box.append([x, y])  # triangle1.3 bas droite
        x = x3 - bl4 - bl3
        x21 = x
        y = y1
        box.append([x, y])  # triangle2.1 bas gauche
        xptm2 = x21 - distx2m
        yptm2 = y1 - disty2m
        x = xptm2
        y = yptm2
        box.append([x, y])  # triangle2m bas gauche

        x = x3 - bl4 - bl3 - bl2
        y = y1
        box.append([x, y])  # p4 bas gauche

        # y = y + bb
        # box.append([x, y])  # pa haut gauche
        # x = x + bl2
        # box.append([x, y])  # ph
        # y = y -bb
        # box.append([x, y])  # pg
        # x = x + bl3
        # box.append([x, y])  # pj
        # y = y + bb
        # box.append([x, y])  # pi
        # x = x + bl4
        # box.append([x, y])  # pb haut gauche
        # y = y - bb
        # box.append([x, y])  # p2 haut droite
        # y = - d / 2
        # box.append([x, y])  # p3 bas droite
        # box.append([x, y])  # pc bas droite
        # x = x- bl4
        # box.append([x, y])  # pu
        # y = y +bb
        # box.append([x, y])  # pv
        # x = x - bl3
        # box.append([x, y])  # px
        # y = y - bb
        # box.append([x, y])  # pw
        # x = x - bl2
        # box.append([x, y])  # pd bas droite
        # y = y + bb
        # box.append([x,y])  # p4 bas gauche
        # bout bas gauche
        for i in range(179):
            x = -a1 - rayon2 * np.sin(np.radians(i))
            y = -a2 + rayon2 * np.cos(np.radians(i))
            box.append([x, y])
        # cutout bas gauche
        for i in range(11, 36):
            x = -dx + rayon * np.sin(np.radians(i)) - dc
            y = + rayon * np.cos(np.radians(i)) - 3.6 * dy
            box.append([x, y])
        x = -(dx - 180)
        y = (-dy)
        box.append([x, y])
        # for i in range(b, theta + 1):
        #  x = -dx - dy * np.sin(np.radians(i))
        #  y = dy * np.cos(np.radians(i))
        #  box.append([x, y])

        # x = -10
        # y = -10
        # box.append([x,y])
        # x = 10
        # box.append([x, y])
        # y = 10
        # box.append([x, y])
        # x = -10
        # box.append([x, y])

    return box



