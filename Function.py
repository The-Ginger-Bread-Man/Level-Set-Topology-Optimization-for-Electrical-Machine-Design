# import files
import os
import shutil
from shapely import geometry
import numpy as np
from numba import njit
from matplotlib import pyplot as plt
from matplotlib import transforms as tran
import scipy


# Internal Function looks for features within other features
# Input level sets
# Outputs list of LS that are within other LS in the format of [exterior LS,interior LS]
def internal(shapes):
    x = []
    y = []
    # Collects points from the LS curves
    for index, element in enumerate(shapes):
        x.append([])
        y.append([])
        points = shapes[index]
        for points_index, points_element in enumerate(points):
            x[index].append(points_element[0])
            y[index].append(points_element[1])
    #print(x)
    container = []
    store = 0
    # compares LS curves together looking for LS's in LS's
    for index, element in enumerate(x):
        # plt.plot(X[a], Y[a],"-")
        if len(x[index]) < 4:  # removes very small LS curves
            container.append([index])
            store += 1
        else:
            n = 0
            for index_2, element_2 in enumerate(x):
                if len(x[index]) < 4:  # removes very small LS curves
                    pass
                else:
                    bool_operator = False
                    if index_2 != index:  # for LS that are not the same
                        poly1 = geometry.Polygon(zip(x[index], y[index]))  # Surface
                        poly2 = geometry.Polygon(zip(x[index_2], y[index_2]))  # Feature
                        bool_operator = poly1.contains(
                            poly2
                        )  # checks surface contains feature (True/False)
                    if bool_operator:
                        if n == 0:
                            # plt.plot(X[b], Y[b],"o")
                            container.append([index, index_2])
                            n += 1
                            store += 1
                        else:
                            # plt.plot(X[b], Y[b],"o")
                            container[store - 1].append(index_2)
                            n += 1
    # plt.show()
    if len(container) > 0:
        container = sorted(container[:], key=len)
        sorting = True
        n = 0
        sorted_list = []
        while sorting:
            test = container.pop(0)
            count = 0
            for element in container:
                for value in element:
                    if test[0] == value:
                        count += 1
            sorted_list.append([count, test])
            n += 1
            if len(container) == 0:
                sorting = False

        sorted_list = sorted(
            sorted_list[:], key=lambda sorted_list: sorted_list[0], reverse=True
        )

        ammended = sorted_list
        for index, element in enumerate(sorted_list):
            for value in range(len(element[1]) - 1):
                test_case = element[1][value + 1]
                for index_2, element_2 in enumerate(ammended):
                    remove = []
                    if index_2 != index:
                        for value_2 in range(len(element_2[1]) - 1):
                            if element_2[1][value_2 + 1] == test_case:
                                remove.append(value_2 + 1)
                        if len(remove) > 0:
                            for value_2 in range(len(remove)):
                                ammended[index_2][1].pop(remove[value_2] - value_2)
        return ammended
    else:
        return []


# Generates random holes in rotor
@njit
def generate_holes(holes, rotor_radius, shaft_radius, angle):
    centers = np.zeros(shape=(holes, 2))
    for index in range(holes):
        center_radius = (
            np.random.random() * (rotor_radius - shaft_radius) + shaft_radius
        )
        center_angle = np.random.random() * (angle * 2)
        centers[index][0] = center_radius
        centers[index][1] = center_angle
    return centers


@njit
def speed_up(
    x,
    y,
    wave_rad,
    mesh_values,
    point_value,
    maximum_b0,
    maximum_b90,
    maximum_stress,
    prev,
    torque,
    best_torque_value,
    best_torque,
    run_torque,
    radius,
    ripple_value,
):
    delta = 2
    for index_mv, value_mv in enumerate(mesh_values):
        for index_x, value_x in enumerate(x):
            for index_y, value_y in enumerate(value_x):
                ratio = 0.5 + (
                    0.5
                    * (((x[index_x][index_y]) ** 2 + (y[index_x][index_y]) ** 2) ** 0.5)
                    / radius
                )
                # print(ratio)
                x_down = x[index_x][index_y] - (wave_rad * ratio)
                x_up = x[index_x][index_y] + (wave_rad * ratio)
                y_down = y[index_x][index_y] - (wave_rad * ratio)
                y_up = y[index_x][index_y] + (wave_rad * ratio)
                if x_down < value_mv[0] < x_up and y_down < value_mv[1] < y_up:
                    point_value[index_x][index_y] += 0.9*((value_mv[2] / maximum_b0) / (
                        value_mv[3] / maximum_b90) - 1) + 0.1*(((-value_mv[4]+maximum_stress )/ maximum_stress)) # 0.8,0.2
                    ripple_value[index_x][index_y] += 0
                    if value_mv[4] > 1e7:
                        point_value[index_x][index_y] = 0
                if torque**2 < best_torque_value**2 and run_torque:
                    point_value[index_x][index_y] += best_torque[index_x][index_y] * (torque**2/best_torque_value**2)**0.5
                if index_mv == len(mesh_values) - 1:
                    point_value[index_x][index_y] += delta
                    if point_value[index_x][index_y] < 0:
                        point_value[index_x][index_y] = 0
                    elif point_value[index_x][index_y] > 2 * delta:
                        point_value[index_x][index_y] = 2 * delta
    for spin in range(1):
        point_value = LinearRegression(
            prev[:], point_value[:], len(x), len(x[0]), 2 * delta
        )
    for index_x, value_x in enumerate(x):
        for index_y, value_y in enumerate(value_x):
            if point_value[index_x][index_y] < 0:
                point_value[index_x][index_y] = 0
            elif point_value[index_x][index_y] > 2 * delta:
                point_value[index_x][index_y] = 2 * delta
            point_value[index_x][index_y] = point_value[index_x][index_y] - delta
    if torque**2 < best_torque_value**2:  # and run_torque:
        return point_value, best_torque_value,ripple_value  # point_value_prev, best_torque_value
    else:
        return point_value, torque,ripple_value
    # return point_value, torque

# Generating initial funtion space
@njit
def xy_coordinates(
    radius_range,
    angle_range,
    angle,
    centers,
    hole_radius,
    list_creator,
    holes,
    shaft_radius,
    rotor_radius,
    delta,
    hole_variaiton = 0.4,
):
    sample = np.zeros(shape=(radius_range * angle_range, 4))
    x = np.copy(list_creator)
    y = np.copy(list_creator)
    z_cont = np.copy(list_creator)
    #Calculating X,Y position in funciton space
    for index_r in range(radius_range):
        for index_a in range(angle_range):
            position = index_r * angle_range + index_a
            sample[position] = [
                shaft_radius
                + ((rotor_radius - shaft_radius) / radius_range) * (index_r + 0.5),
                (-angle) + ((2 * angle) / angle_range) * (index_a + 0.5),
                0,
                0,
            ]
            x[index_r][index_a] = sample[position][0] * np.sin(sample[position][1])
            y[index_r][index_a] = sample[position][0] * np.cos(sample[position][1])
    # adds holes to fucntion space based on distance
    for index, element in enumerate(sample):
        sign = 1
        size = 4
        for value in range(holes):
            x_test = element[0] * np.sin(element[1])
            y_test = element[0] * np.cos(element[1])
            c_x = centers[value][0] * np.sin(centers[value][1] - angle)
            c_y = centers[value][0] * np.cos(centers[value][1] - angle)
            hole_radius_modified = hole_radius * (1-hole_variaiton) + hole_variaiton * hole_radius * (
                (centers[value][0] - shaft_radius) / (rotor_radius - shaft_radius)
            ) # varies hole radius with radial position 40% variation
            # distance from boundary
            if (
                x_test**2 + y_test**2 - 2 * (c_x * x_test + c_y * y_test)
                < hole_radius_modified**2 - c_x**2 - c_y**2
            ):
                sign = 0
                size = -(hole_radius_modified**2) + (
                    x_test**2
                    + y_test**2
                    - 2 * (c_x * x_test + c_y * y_test)
                    + c_x**2
                    + c_y**2
                )  # * 1e6
            if size > -(hole_radius_modified**2) + (
                x_test**2
                + y_test**2
                - 2 * (c_x * x_test + c_y * y_test)
                + c_x**2
                + c_y**2
            ):
                size = -(hole_radius_modified**2) + (
                    x_test**2
                    + y_test**2
                    - 2 * (c_x * x_test + c_y * y_test)
                    + c_x**2
                    + c_y**2
                )  # * 1e6
        # sets limit from boundary
        if size > delta:
            size = delta
        elif size < -delta:
            size = -delta
        sample[index][2] = sign
        sample[index][3] = size
    # Generates full function space
    for index_r in range(radius_range):
        for index_a in range(angle_range):
            if (
                index_r == 0
                or index_r == radius_range - 1
            ):
                z_cont[index_r][index_a] = delta
            elif (
                index_a == 0 or
                index_a == angle_range - 1
            ):
                z_cont[index_r][index_a] = delta
            else:
                z_cont[index_r][index_a] = sample[index_r * angle_range + index_a][3]
    return z_cont, x, y


def get_data(length, angle1_lines, angle2_lines,stress_lines):
    mesh_values = np.zeros(shape=(length, 5))
    maximum_stress = 1
    maximum_b0 = 1
    maximum_b90 = 1
    print("\rProcessing Data", end="")
    for index in range(length):
        try:
            results = str(angle1_lines[index]).split("(")
            results90 = str(angle2_lines[index]).split("(")
            Resultsstress = str(stress_lines[index]).split("(")
            if len(results) == 2:
                results = results[1].split(")")
                results90 = results90[1].split(")")
                Resultsstress = Resultsstress[1].split(")")
                positions = results[0].split(",")
                value = float(results[1].split(",")[1])
                value90 = float(results90[1].split(",")[1])
                vlaueStress = float(Resultsstress[1].split(",")[1])
                xav = (
                    (float(positions[0]) + float(positions[3]) + float(positions[6]))
                    * 10**3
                    / 3
                )
                yav = (
                    (float(positions[1]) + float(positions[4]) + float(positions[7]))
                    * 10**3
                    / 3
                )
                mesh_values[index][0] = xav
                mesh_values[index][1] = yav
                mesh_values[index][2] = value
                mesh_values[index][3] = value90
                mesh_values[index][4] = vlaueStress  # Value stress
                if maximum_b0 < value:
                    maximum_b0 = value
                if maximum_b90 < value90:
                    maximum_b90 = value90
                if maximum_stress < vlaueStress:
                    maximum_stress = vlaueStress
        except:
            pass
    return mesh_values, maximum_stress, maximum_b0, maximum_b90


def get_file(directory, angle1, angle2):
    angle1Open = open(directory + "b" + str(angle1) + ".pos")
    angle1_lines = angle1Open.readlines()
    angle1Open.close()
    angle2Opne = open(directory + "b" + str(angle2) + ".pos")
    angle2_lines = angle2Opne.readlines()
    angle2Opne.close()
    torque_file = open(directory + "T" + str(angle1) + ".pos")
    torque_value = float(torque_file.readlines()[0].split(" ")[2].split("\\")[0])
    stressOpen = open(directory + "SigMis.pos")
    stress_lines = stressOpen.readlines()
    stressOpen.close()
    return angle1_lines, angle2_lines, torque_value,stress_lines

def get_ripple_file(directory, loops):
    T_Ripple = []
    for value in range(loops):
        file_read = open(directory + "T_rip_.pos" + str(value))
        T_Ripple.append(file_read.readlines())
        file_read.close()
    return T_Ripple
def get_ripple_data(length, files,multiply):
    mesh_values = np.zeros(shape=(length, 3))
    range_data = np.zeros(shape=(len(files), 1))
    print("\rProcessing Data", end="")
    for index in range(length):
        try:
            B_max = 0
            B_min = 10
            for a in range(len(files)):
                B_test = float(str(files[a][index]).split("(")[1].split(")")[1].split(",")[1])
                range_data[a] += B_test
                if B_max < B_test:
                    B_max = B_test
                elif B_min > B_test:
                    B_min = B_test
            results = str(files[0][index]).split("(")[1].split(")")
            positions = results[0].split(",")
            value = (B_max - B_min)*multiply
            xav = (
                (float(positions[0]) + float(positions[3]) + float(positions[6]))
                * 10**3
                / 3
            )
            yav = (
                (float(positions[1]) + float(positions[4]) + float(positions[7]))
                * 10**3
                / 3
            )
            mesh_values[index][0] = xav
            mesh_values[index][1] = yav
            mesh_values[index][2] = value
        except:
            pass
    return mesh_values, range_data

#draws rotor and shaft holes
def circle_rotor_shaft():
    rad = plt.Circle((0, 0), 49.5, color="black", fill=False)
    shaft = plt.Circle((0, 0), 12.5, color="black", fill=False)
    plt.gca().add_patch(rad)
    plt.gca().add_patch(shaft)

#Rotates images
def rotate():
    base = plt.gca().transData
    rot_0 = tran.Affine2D().rotate_deg(0)
    rot_90 = tran.Affine2D().rotate_deg(90)
    rot_180 = tran.Affine2D().rotate_deg(180)
    rot_270 = tran.Affine2D().rotate_deg(270)
    return base, rot_0, rot_90, rot_180, rot_270

# draws funciotn space
def draw_curve(x, y, z_cont, min_length, x1, y1, iteration, gausian_filter,File):
    ls_contours = []
    contours = 0
    #draws fucntion space at 0 level
    w = scipy.signal.windows.gaussian(51, gausian_filter)
    z_cont = scipy.signal.sepfir2d(z_cont, w, w)
    base, rot_0, rot_90, rot_180, rot_270 = rotate()
    ls = plt.contour(x, y, z_cont, [0], transform=rot_0 + base)
    ls_90 = plt.contour(x, y, z_cont, [0], transform=rot_90 + base)
    ls_180 = plt.contour(x, y, z_cont, [0], transform=rot_180 + base)
    ls_270 = plt.contour(x, y, z_cont, [0], transform=rot_270 + base)
    circle_rotor_shaft() # draws rotor and shaft
    # collects 0 level set
    ls_contours_test = ls.allsegs[0]
    plt.close()
    # Presents fucniton space to user
    base, rot_0, rot_90, rot_180, rot_270 = rotate()
    for element in ls_contours_test:
        v = element
        if len(v[:, 0]) > min_length:
            ls_contours.append(element)
            contours += 1
            plt.plot(v[:, 0], v[:, 1], "r", transform=rot_0 + base)
            plt.plot(v[:, 0], v[:, 1], "r", transform=rot_90 + base)
            plt.plot(v[:, 0], v[:, 1], "r", transform=rot_180 + base)
            plt.plot(v[:, 0], v[:, 1], "r", transform=rot_270 + base)
    circle_rotor_shaft()
    plt.savefig(File+"Images\\loopnumber" + str(iteration) + ".pdf")
    plt.savefig(File+"Images\\loopnumber" + str(iteration) + ".png")
    plt.close()
    return ls_contours, contours

@njit
def cleanup(z_cont, radius_range, angle_range, delta):
    for index_r in range(radius_range):
        for index_a in range(angle_range):
            if (
                index_r == 0
                or index_r == radius_range - 1
            ):
                z_cont[index_r][index_a] = delta
            elif (
                index_a == 0 or
                index_a == angle_range - 1
            ):
                z_cont[index_r][index_a] = delta
            elif z_cont[index_r][index_a] > delta:
                z_cont[index_r][index_a] = delta
            elif z_cont[index_r][index_a] < -delta:
                z_cont[index_r][index_a] = -delta
    return z_cont

#@njit
def update_2(point_value, dt, lag,GRAD, z_cont, radius_range, angle_range, delta):
    v = np.zeros(shape=(radius_range, angle_range))
    for index_r in range(radius_range):
        for index_a in range(angle_range):
            if point_value[index_r][index_a] != 0 and (point_value[index_r][index_a])**2 != delta**2:
                #print("working",point_value[index_r][index_a])
                v[index_r][index_a] = (point_value[index_r][index_a]) * dt - GRAD*(lag[index_r][index_a]-delta)  # Lag
            else:
                v[index_r][index_a] = 0

    #z_cont = np.subtract(z_cont, v)
    for index_r in range(radius_range):
        for index_a in range(angle_range):
            z_cont[index_r][index_a] = z_cont[index_r][index_a] - v[index_r][index_a]
    z_cont = cleanup(z_cont, radius_range, angle_range, delta)

    return z_cont

def update(point_value, dt, lag, z_cont, radius_range, angle_range, delta):
    v = np.zeros(shape=(radius_range, angle_range))
    for index_r in range(radius_range):
        for index_a in range(angle_range):
            if point_value[index_r][index_a] != 0 and (point_value[index_r][index_a])**2 != delta**2:
                #print("working",point_value[index_r][index_a])
                v[index_r][index_a] = (point_value[index_r][index_a]) * dt - lag[index_r][index_a]  # Lag
            else:
                v[index_r][index_a] = 0

    #z_cont = np.subtract(z_cont, v)
    for index_r in range(radius_range):
        for index_a in range(angle_range):
            z_cont[index_r][index_a] = z_cont[index_r][index_a] - v[index_r][index_a]
    z_cont = cleanup(z_cont, radius_range, angle_range, delta)

    return z_cont


@njit
def volume_calc(z_cont, radius_range, angle_range):
    volume = 0
    for index_r in range(radius_range):
        for index_a in range(angle_range):
            if z_cont[index_r][index_a] > 0:
                volume += 1
    return volume

#Shrink mesh
def minimise_msh(directory, file, mag, type):
    file_read = open(str(directory) + str(file) + str(type)).readlines()
    file_new = open(str(directory) + str(file) + "_mini" + str(type), "w")
    if type == ".msh":
        for index, element in enumerate(file_read):
            test_case = file_read[index].split(" ")
            if not len(test_case) > 1 or index < 3:
                file_new.write(file_read[index])
            elif float(test_case[1]) != int(float(test_case[1])) or float(
                    test_case[0]
            ) != int(float(test_case[0])):
                modified = []
                for index_tc, element_tc in enumerate(test_case):
                    try:
                        if len(test_case[index_tc].split(".")) > 1:
                            if float(element_tc) != int(float(element_tc)):
                                modified.append(str(float(element_tc) * 10 ** mag))
                            else:
                                modified.append(str(float(element_tc) * 10 ** mag))
                        else:
                            if (
                                    (element_tc == "100" or element_tc == "50")
                                    and 0 < index_tc < 10
                                    and len(test_case) < 10
                            ):
                                modified.append(str(float(element_tc) * 10 ** mag))
                            else:
                                modified.append(str(element_tc))
                    except:
                        if (str(element_tc) == "100" or str(element_tc) == "50") and len(test_case) < 10:
                            modified.append(str(float(element_tc) * 10 ** mag))
                        else:
                            modified.append(element_tc)
                modified = " ".join(modified)
                file_new.write(modified)
            else:
                file_new.write(file_read[index])
    else:
        for index, element in enumerate(file_read):
            try:
                test_case = file_read[index].split("(")
                mesh_element = test_case[1].split(")")
                coordinates = mesh_element[0].split(",")
                new_coordinates = "("
                for position, coordinate in enumerate(coordinates):
                    if position == len(coordinates) - 1:
                        new_coordinates += str(float(coordinate) * 10**mag) + str(")")
                    else:
                        new_coordinates += str(float(coordinate) * 10**mag) + str(",")
                file_new.write(test_case[0])
                file_new.write(new_coordinates)
                file_new.write(mesh_element[1])
            except:
                file_new.write(element)
    file_new.close()

#removes some points
def filter(level_sets, density):
    new_level_sets = []
    for index, element in enumerate(level_sets):
        new_level_sets.append(element[::density])
    return new_level_sets

def printpoint(b):
    obj = dict(zip(b[0::2], b[1::2]))
    #print(obj)
    if obj['0'] == 'LINE':
        #print('{},{},{},{}'.format(obj['10'], obj['20'], obj['11'], obj['21']))
        return [float(obj['10']), float(obj['20']), float(obj['11']), float(obj['21'])]
def DXF_import(file):
    input = open(file).readlines()
    buffer = ['0', 'fake']    # give first pass through loop something to process
    line_collection = []
    loop_collection = []
    for line in input:
        line = line.rstrip()
        if line == '0':         # we've started a new section, so
            line_test=printpoint(buffer)      # handle the captured section
            print(line_test)
            if line_test != None:
                plt.plot(line_test[0],line_test[1],"o")
                if len(line_collection) < 1:
                    line_collection.append(line_test)
                elif round(float(line_collection[len(line_collection)-1][2]),2) == round(float(line_test[0]),2):
                    line_collection.append(line_test)
                    plt.plot(line_test[0],line_test[1],"x")
                    #print(line_test)
                elif  round(float(line_collection[len(line_collection)-1][0]),2) == round(float(line_test[2]),2):
                    #print(line_test,"1")
                    plt.plot(line_test[0],line_test[1],"x")
                    #line_collection.append([line_test[2],line_test[3],line_test[0],line_test[1]])
                    line_collection.append(line_test)
                elif  round(float(line_collection[len(line_collection)-1][0]),2) == round(float(line_test[0]),2):
                    line_collection.append([line_test[2],line_test[3],line_test[0],line_test[1]])
                    plt.plot(line_test[2],line_test[3],"x")
                    #print(line_test,"2")
                    line_collection.append(line_test)
                elif  round(float(line_collection[len(line_collection)-1][2]),2) == round(float(line_test[2]),2):
                    line_collection.append([line_test[2],line_test[3],line_test[0],line_test[1]])
                    plt.plot(line_test[2],line_test[3],"x")
                    #print(line_test,"3")
                else:
                    #print(line_test)
                    #print(line_collection)
                    line_collection.append([line_collection[-1][2],line_collection[-1][3],line_collection[0][0],line_collection[0][1]])
                    loop_collection.append(line_collection)
                    line_collection = [line_test]
            buffer = []             # and start a new one
        buffer.append(line)
    line_collection.append([line_collection[-1][2],line_collection[-1][3],line_collection[0][0],line_collection[0][1]])
    loop_collection.append(line_collection)
    #remove rotor boundaries
    #print(loop_collection)
    loops = []
    angle = 45*np.pi/180
    #rotate_anlges=[np.cos(angle),-np.sin(angle)]
    #Rotor _ opitimized
    #loop_collection.pop(len(loop_collection)-1)
    #loop_collection.pop(len(loop_collection)-1)
    #loop_collection.pop(len(loop_collection)-1)
    # conventional
    #loop_collection.pop(0)
    #loop_collection.pop(0)
    #loop_collection.pop(0)

    for loop in range(len(loop_collection)):#-3
        loops.append([])
        for point in loop_collection[loop]:
            loops[loop].append([point[2]*np.cos(angle)+point[3]*-np.sin(angle),
                                point[2]*np.sin(angle)+point[3]*np.cos(angle)])
            plt.plot(point[2],point[3],"x")
    plt.savefig("points.pdf")
    #plt.show()
    return loops

@njit
def LinearRegression(A,B,rows,columbs,delta):
    Atemp = np.zeros((rows, columbs))
    for r in range(rows):
        for c in range(columbs):
            P1 = A[EdgeFinder(r , rows)][EdgeFinder(c - 1, columbs)]
            P2 = A[EdgeFinder(r , rows)][EdgeFinder(c + 1, columbs)]
            P3 = A[EdgeFinder(r + 1, rows)][EdgeFinder(c , columbs)]
            P4 = A[EdgeFinder(r - 1, rows)][EdgeFinder(c , columbs)]
            Magnitudes = np.array([B[EdgeFinder(r , rows)][EdgeFinder(c - 1, columbs)],
                          B[EdgeFinder(r , rows)][EdgeFinder(c + 1, columbs)],
                          B[EdgeFinder(r + 1, rows)][EdgeFinder(c , columbs)],
                          B[EdgeFinder(r - 1, rows)][EdgeFinder(c , columbs)]])
            if P1 > 0 and P2 > 0 and P3 > 0 and P4 > 0:
                Atemp += 0
            elif P1 < 0 and P2 <0 and P3 < 0 and P4 < 0:
                Atemp += 0
            else:
                P1,P2,P3,P4 = Update(P1,P2,P3,P4,Magnitudes)
                Atemp[r][c] += P1 * 0.25 + P2 * 0.25 + P3 * 0.25 + P4 * 0.25
    return Atemp
@njit
def EdgeFinder(position,limits):
    if position < 0:
        return limits-1
    elif position == limits:
        return 0
    else:
        return position
@njit
def  Update(P1,P2,P3,P4,Mags):
    Points = np.array([P1,P2,P3,P4])
    xy = np.array([[0,0],[1,0],[0,1],[1,1]])
    if P1 > 0:
        if P2 > 0:
            if P3 > 0: # p1,p2,p3 + p4 -
                Corners = np.array([[1,3,2,3]])
            else:
                if P4 > 0: # p1,p2,p4 + p3 -
                    Corners = np.array([[0,2,3,2]])
                else:# p1,p2 + p3,p4  -
                    Corners = np.array([[0,2,1,3]])
        else:
            if P3> 0:
                if P4> 0:# p3,p4,p1 + p2  -
                    Corners = np.array([[0,1,3,1]])
                else:# p3,p1 + p4,p2  -
                    Corners = np.array([[0,1,2,3]])
            else:
                if P4> 0:# p1,p4 + p3,p2  -
                    Corners = np.array([[0,1,0,2],[3,1,3,2]])
                else:# p1 + p3,p4,p2  -
                    Corners = np.array([[0,2,0,2]])
    else:
        if P2 > 0:
            if P3> 0:
                if P4> 0: # p4,p2,p3 + p1 -
                    Corners = np.array([[1,0,2,0]])
                else: # p2,p3 + p4,p1 -
                    Corners = np.array([[1,0,2,0],[1,3,2,3]])
            else:
                if P4> 0: # p4,p2 + p1,p3 -
                    Corners = np.array([[1,0,2,0]])
                else: # p2 + p4,p1,p3 -
                    Corners = np.array([[1,0,1,3]])
        else:
            if P3> 0:
                if P4> 0: # p4,p3 + p1,p2 -
                    Corners = np.array([[2,0,3,1]])
                else: # p3 + p4,p2,p1 -
                    Corners = np.array([[1,0,2,0]])
            else: # p4 + p1,p2,p3 -
                Corners = np.array([[3,1,3,2]])
    Shift=0
    m=0
    for a in Corners:
        Shift += (Mags[a[0]] + Mags[a[1]]+ Mags[a[2]]+ Mags[a[3]])/4
        Shift = Shift/(m+1)
        m+=1
    #if type(Corners[1]) != int:
    #    Shift = (Mags[Corners[0][0]] + Mags[Corners[0][1]]+ Mags[Corners[0][2]]+ Mags[Corners[0][3]])/4
    #    Shift += (Mags[Corners[1][0]] + Mags[Corners[1][1]]+ Mags[Corners[1][2]]+ Mags[Corners[1][3]])/4
    #    return P1+Shift,P2+Shift,P3+Shift,P4+Shift
    #else:
    #    Interp_A = Points[Corners[0]]/(Points[Corners[0]]-Points[Corners[1]])
    #    Interp_B = Points[Corners[2]]/(Points[Corners[2]]-Points[Corners[3]])
    #    Position_A = [Interp_A*xy[Corners[0]][0],Interp_A*xy[Corners[0]][1]]
    #    Position_B = [Interp_B*xy[Corners[2]][0],Interp_B*xy[Corners[2]][1]]
    #    Grad = (Position_A[1]-Position_B[1])/(Position_A[0]-Position_B[0])
    #    Shift = (Mags[Corners[0]] + Mags[Corners[1]]+ Mags[Corners[2]]+ Mags[Corners[3]])/4
    #    return P1+Shift,P2+Shift,P3+Shift,P4+Shift
    return P1+Shift,P2+Shift,P3+Shift,P4+Shift


def copy_file(source_path, destination_path):
    """
    Copies a file from source_path to destination_path.
    Overwrites the destination if it already exists.
    """
    try:
        # Validate that the source file exists
        if not os.path.isfile(source_path):
            raise FileNotFoundError(f"Source file not found: {source_path}")

        # Ensure the destination directory exists
        dest_dir = os.path.dirname(destination_path)
        if dest_dir and not os.path.exists(dest_dir):
            os.makedirs(dest_dir)  # Create directories if needed

        # Perform the copy
        shutil.copy2(source_path, destination_path)  # copy2 preserves metadata
        #print(f"File copied successfully from '{source_path}' to '{destination_path}'.")

    except FileNotFoundError as fnf_error:
        print(f"Error: {fnf_error}")
    except PermissionError:
        print("Error: Permission denied. Check file and folder permissions.")
    except IsADirectoryError:
        print("Error: Destination path is a directory, not a file.")
    except Exception as e:
        print(f"Unexpected error: {e}")

def make_run(Files):
    for a in range(Files):
        try:
            os.mkdir("C:\\Users\\fh18930\\PycharmProjects\\Simple_Example\\Run_Multi\\" + str(a))
            os.mkdir("C:\\Users\\fh18930\\PycharmProjects\\Simple_Example\\Run_Multi\\" + str(a) + "\\images")
            os.remove(
                "C:\\Users\\fh18930\\PycharmProjects\\Simple_Example\\Run_Multi\\" + str(a) + "\\LSModel_mini.pro")
        except:
            pass
        try:
            if a > 0:
                copy_file(
                    "C:\\Users\\fh18930\\PycharmProjects\\Simple_Example\\Run_Multi\\0\\LSModel_mini.pro",
                    "C:\\Users\\fh18930\\PycharmProjects\\Simple_Example\\Run_Multi\\" + str(a) + "\\")
        except:
            pass