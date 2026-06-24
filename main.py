##############Set UP#################
import time
from logging import error
from os import abort
import numpy as np
import Model
import Function

print("")
print("")
print("Begining Topology Optimization")
print("")
print("")

# File Directory
directory = "C:\\Users\\fh18930\\PycharmProjects\\Optimization_Publication\\Run\\"
GetDP = "C:\\Users\\fh18930\\OneDrive - University of Bristol\\Desktop\\onelab-Windows64\\getDP.exe"

# Variables for conventional optimizaiton
constant_start = False
save_start = False
max_generation = 100            # Generation Maximum
gausian_filter = 0.25           # Gaussian Filter
radius_range = 50               # Points radialy in funciton space
angle_range = radius_range      # Points angularly in funciton space
dt = 0.2                        # Time step reduce [s]
delta = 2                       # Controlls bounds of LS funciton
wave_rad = 1                    # Radius Filter [mm]
holes = 500                     # Initial holes in model
hole_radius = 1                 # Initial hole radius [mm]

# Machine Dimensions
current = 2.58                  # Current Load [A]
rotor_radius = 49.5             # Rotor radius [mm]
shaft_radius = 15               # Shaft radius [mm]
stator_radius = 80
tooth_width = 6
airgap = 0.5
shaft_buffer = 2.5
poles = 4
airgap_mesh = 0.1
rotor_mesh = 1
back_iron = 13
stator_mesh = 5

# Variables for torque ripple constraint
run_ripple = False
wind_up = 0

# Advanced variables
skip = 1
angle1 = 0
angle2 = 45 * np.pi / 180
angle = np.pi / 4  # Rotor angle modeled
min_length = 8
best_torque_ripple_value = 100
torque_ripple = np.zeros(shape =(1,3))
torque_ripple[0] = 100
dryrun = False
distance_calculation = False

############################################
############## Setup #######################
############################################
processing = True
sim_timer = time.time()             # Starts timer
iteration = 1
repeat = 0
lag = 0
maximum_rip = 1000
historical_volume = []

if constant_start and save_start:
    error("please set one to constant_start or save_start to False")
    abort()

###############Level Set Variables##################

# draws surfaces of holes
centers = Function.generate_holes(holes, rotor_radius, shaft_radius, angle)
# Draws surface of rotor
x1 = []
y1 = []
for a in range(51):
    x1.append(rotor_radius * np.sin((25 - a) * (angle / 25)))
    y1.append(rotor_radius * np.cos((25 - a) * (angle / 25)))
for a in range(51):
    x1.append(shaft_radius * np.sin((-25 + a) * (angle / 25)))
    y1.append(shaft_radius * np.cos((-25 + a) * (angle / 25)))
x1.append(rotor_radius * np.sin(angle))
y1.append(rotor_radius * np.cos(angle))

# Generates funciton space
listCreator = np.zeros(shape=(radius_range, angle_range))
z_cont, x, y = Function.xy_coordinates(
    radius_range,
    angle_range,
    angle,
    centers[:],
    hole_radius,
    listCreator,
    holes,
    shaft_radius,
    rotor_radius,
    delta,
)

# Saves initial Condition
if save_start:
    constant_start = open("Constant_start.pos", "w")
    for a in range(len(z_cont)):
        for b in range(len(z_cont[a])):
            if b == len(z_cont[a]) - 1:
                constant_start.write(str(z_cont[a][b]) + ",\n")
            else:
                constant_start.write(str(z_cont[a][b]) + ",")
    constant_start.close()

# Loads initial Condition
if constant_start:
    constant_start = open("Constant_start.pos", "r").readlines()
    for a in range(len(constant_start)):
        list = constant_start[a].split(",")
        for b in range(len(list) - 1):  # removes \n
            z_cont[a][b] = float(list[b])

z_cont_new = listCreator.copy()
z_cont_old = z_cont_new.copy()
z_cont_new = z_cont.copy()

# Draw Curve
ls_contours, contours = Function.draw_curve(
    x, y, z_cont, min_length, x1, y1, 0, gausian_filter,directory
)

# Reduce number of points
ls_contours = Function.filter(ls_contours, skip)

# builds models
container_2 = Function.internal(ls_contours)
if not dryrun:
    Model.Model(
        airgap,
        rotor_radius,
        stator_radius,
        tooth_width,
        back_iron,
        shaft_radius-shaft_buffer,
        angle,
        airgap_mesh,
        rotor_mesh,
        stator_mesh,
        angle1,
        angle2,
        ls_contours,
        container_2,
        directory,
        current,
        GetDP,
    )

# Collecting Data
point_value_setup = np.zeros(shape=(len(x), len(x[0])))
best_torque_value = 0
oldtorque = 0
best_torque = point_value_setup
# Get results
angle1_lines, angle2_lines, torque,stress_lines = Function.get_file(directory, angle1, angle2)
newtorque = torque
length = len(angle1_lines)
mesh_values, maximum_stress, maximum_b0, maximum_b90 = Function.get_data(
    length, angle1_lines, angle2_lines,stress_lines
)
# Carry out data processing
new_point_value = best_torque.copy()
if run_ripple:
    T_ripple = Function.get_ripple_file(directory,4)
    torque_ripple = np.zeros(shape=[4, 1])
    for a in range(len(T_ripple)):
        torque_ripple[a] = T_ripple[a][0].split(" ")[-1].split("\n")[0]

point_value, best_torque_value,lag = Function.speed_up(
    x,
    y,
    wave_rad,
    mesh_values,
    np.copy(point_value_setup),
    maximum_b0,
    maximum_b90,
    maximum_stress,
    np.copy(z_cont),
    torque,
    best_torque_value,
    best_torque,
    False,
    rotor_radius,
    np.copy(point_value_setup),
)

old_point_value = new_point_value.copy()
new_point_value = point_value.copy()
range_energy = 1000
best_torque_ripple = z_cont.copy()
while processing:
    #gausian_filter = 0.1 + iteration/(50)
    print("\rRunning iteration:" + str(iteration), end="")
    print(" ")
    historical_volume.append(Function.volume_calc(z_cont, radius_range, angle_range))

    if iteration > max_generation:
        processing = False

    #print(lag, (historical_volume[iteration - 1]) / (radius_range * angle_range))
    print("\rCurrent Torque : " + str(round(torque,1))+ ", best Torque : " + str(round(best_torque_value,1)), end="")
    print(" ")
    old_range_energy = range_energy
    range_energy = ((torque_ripple.max()-torque_ripple.min())/torque_ripple.mean())#maximum_rip.max() - maximum_rip.min()
    if range_energy < best_torque_ripple_value and iteration > 5:
        best_torque_ripple_value = range_energy
        best_torque_ripple = z_cont.copy()
        print("\rTorque ripple : " + str(round(range_energy*100,1)), end="")
        print(" ")
    else:
        print("\rTorque ripple : " + str(round(range_energy*100,1))+ " best Torque ripple : "+ str(round(best_torque_ripple_value*100,1)), end="")
        print(" ")

    if torque == best_torque_value:
        best_torque = z_cont.copy()

    if iteration > 10 and run_ripple:
        wind_up += 0.01
        if wind_up > 0.25:
            wind_up = 0.25
    direction = -1
    GRADIENT = wind_up*direction*(best_torque_ripple_value-old_range_energy)/0.1
    z_cont = Function.update_2(
        point_value, dt, best_torque_ripple,GRADIENT, z_cont, radius_range, angle_range, delta
    )
    z_cont_old = z_cont_new.copy()
    z_cont_new = z_cont.copy()
    ls_contours, contours = Function.draw_curve(
        x, y, z_cont, min_length, x1, y1, iteration, gausian_filter,directory
    )

    # Reduce points
    ls_contours = Function.filter(ls_contours, skip)

    container_2 = Function.internal(ls_contours)
    if iteration < max_generation-4:
        distance_calculation = True
    if not dryrun:
        Model.Model(
            airgap,
            rotor_radius,
            stator_radius,
            tooth_width,
            back_iron,
            shaft_radius-shaft_buffer,
            angle,
            airgap_mesh,
            rotor_mesh,
            stator_mesh,
            angle1,
            angle2,
            ls_contours,
            container_2,
            directory,
            current,
            GetDP,
        )

    angle1_lines, angle2_lines, torque,stress_lines = Function.get_file(directory, angle1, angle2)
    length = len(angle1_lines)
    mesh_values, maximum_stress, maximum_b0, maximum_b90 = Function.get_data(
        length, angle1_lines, angle2_lines,stress_lines
    )

    oldtorque = newtorque
    newtorque = torque

    if run_ripple:
        T_ripple = Function.get_ripple_file(directory,4)
        torque_ripple = np.zeros(shape=[4,1])
        for a in range(len(T_ripple)):
            torque_ripple[a] = T_ripple[a][0].split(" ")[-1].split("\n")[0]

    point_value, best_torque_value,lag = Function.speed_up(
        x,
        y,
        wave_rad,
        mesh_values,
        np.copy(point_value_setup),
        maximum_b0,
        maximum_b90,
        maximum_stress,
        np.copy(z_cont),
        torque,
        best_torque_value,
        best_torque,
        False,
        rotor_radius,
        np.copy(point_value_setup),
    )
    old_point_value = new_point_value.copy()
    new_point_value = point_value.copy()

    iteration += 1


Model.Model(
    airgap,
    rotor_radius,
    stator_radius,
    tooth_width,
    back_iron,
    shaft_radius-shaft_buffer,
    angle,
    airgap_mesh,
    rotor_mesh,
    stator_mesh,
    angle1,
    angle2,
    ls_contours,
    container_2,
    directory,
    current,
    GetDP,
)
angle1_lines, angle2_lines, torque, stress_lines = Function.get_file(directory, angle1, angle2)
length = len(angle1_lines)
mesh_values, maximum_stress, maximum_b0, maximum_b90 = Function.get_data(
    length, angle1_lines, angle2_lines, stress_lines
)
oldtorque = newtorque
newtorque = torque
print("\rCurrent Torque : " + str(round(torque,1)) + ", best Torque : " + str(round(best_torque_value,1)), end="")

print("Simulation time of :" + str(sim_timer - time.time()) + "s")
#print(x)
#print(y)
FileLoadFile = open(directory+"Result.txt","w")
FileLoadFile.write("[")
for a in range(len(best_torque)):
    FileLoadFile.write("[")
    for b in range(len(best_torque)):
        FileLoadFile.write(str(best_torque[a][b]) + ",")
    FileLoadFile.write("],\n")
FileLoadFile.write("]")
FileLoadFile.close()
print(z_cont)
