# import files
import time
from logging import error

import gmsh
import numpy as np
import Function


def Model(
        AirGap,
        RotorBore,
        Stator_radius,
        TotalToothWidth,
        BackIron,
        ShaftBore,
        SegmentAngle,
        MeshDensityAG,
        MeshDensityR,
        StatorDensity,
        Angle1,
        Angle2,
        Features,
        Container2,
        directory,
        current,
        getDP,
        Winding_division=8.5,
):
    #start_time = time.time()

    StatorOD = Stator_radius*2
    StatorID = RotorBore * 2 + (2 * AirGap)
    ToothWidth = TotalToothWidth / 2
    spp = 3
    ToothAngleOD2 = np.arcsin((ToothWidth / 2) / ((StatorID) / 2 + 2 + Winding_division)) * 2
    ToothAngleOD = np.arcsin((ToothWidth / 2) / ((StatorID) / 2 + 2)) * 2
    ToothAngleID = np.arcsin((ToothWidth / 2) / ((StatorOD) / 2 - BackIron)) * 2

    slot_opening = 1
    tooth_tip = (StatorOD*np.pi/(spp*4*3)-TotalToothWidth-slot_opening*2)
    #print(tooth_tip)
    #(3.9) * np.pi / 180
    teeth_angle_outer = np.arcsin((ToothWidth/2+tooth_tip) / ((StatorID) + 2))
    teeth_angle_inner = np.arcsin((ToothWidth/2+tooth_tip) / ((StatorID)))
    mag = -3 # building models at m and scale down to mm

    # initialize gmsh
    print("\rgenerating mesh", end="")
    gmsh.initialize()

    # Set up gmsh model
    gmsh.option.setNumber("General.Terminal", 0)                #Turn of terminal messages
    gmsh.option.setNumber("Mesh.SurfaceFaces", 0)               #Surface faces off
    gmsh.option.setNumber("Mesh.SurfaceEdges", 0)               #Surface edges off
    gmsh.option.setNumber("Geometry.Surfaces", 0)               #Surface off
    gmsh.option.setNumber("Geometry.SurfaceLabels", 0)          #Surface labels off
    gmsh.option.setNumber("Geometry.SurfaceType", 0)            #Surface type off
    gmsh.option.setNumber("Geometry.Curves", 1)                 #Curves on
    gmsh.option.setNumber("Geometry.Points", 0)                 #Curves off
    gmsh.model.add("Models\\LSModel")

    # define points
    gmsh.model.geo.addPoint(0, 0, 0, 0, 1)
    gmsh.model.geo.addPoint(
        ShaftBore * np.sin(SegmentAngle),
        ShaftBore * np.cos(SegmentAngle),
        0,
        MeshDensityR,
        2,
    )
    gmsh.model.geo.addPoint(
        -ShaftBore * np.sin(SegmentAngle),
        ShaftBore * np.cos(SegmentAngle),
        0,
        MeshDensityR,
        3,
    )
    gmsh.model.geo.addPoint(
        RotorBore * np.sin(SegmentAngle),
        RotorBore * np.cos(SegmentAngle),
        0,
        MeshDensityAG * 2,
        4,
    )
    gmsh.model.geo.addPoint(
        -RotorBore * np.sin(SegmentAngle),
        RotorBore * np.cos(SegmentAngle),
        0,
        MeshDensityAG * 2,
        5,
    )
    gmsh.model.geo.addPoint(
        (RotorBore * 2 + (AirGap - MeshDensityAG / 2)) / 2 * np.sin(SegmentAngle),
        (RotorBore * 2 + (AirGap - MeshDensityAG / 2)) / 2 * np.cos(SegmentAngle),
        0,
        MeshDensityAG,
        6,
    )
    gmsh.model.geo.addPoint(
        -(RotorBore * 2 + (AirGap - MeshDensityAG / 2)) / 2 * np.sin(SegmentAngle),
        (RotorBore * 2 + (AirGap - MeshDensityAG / 2)) / 2 * np.cos(SegmentAngle),
        0,
        MeshDensityAG,
        7,
    )

    n = 0
    for a in range(3):  # moving boundry points
        n += 1
        gmsh.model.geo.addPoint(
            ((RotorBore * 2 + (AirGap - MeshDensityAG / 2)) / 2)
            * np.sin(SegmentAngle * (3 + (2 * a))),
            (RotorBore * 2 + (AirGap - MeshDensityAG / 2))
            / 2
            * np.cos(SegmentAngle * (3 + (2 * a))),
            0,
            MeshDensityAG,
            11 + n,
        )

    StatorStartPoint = n + 12
    gmsh.model.geo.addPoint(
        StatorOD / 2 * np.sin(np.pi / 4),
        StatorOD / 2 * np.cos(np.pi / 4),
        0,
        StatorDensity,
        StatorStartPoint,
    )
    gmsh.model.geo.addPoint(
        StatorOD / 2 * np.sin(-np.pi / 4),
        StatorOD / 2 * np.cos(-np.pi / 4),
        0,
        StatorDensity,
        StatorStartPoint + 1,
    )
    gmsh.model.geo.addPoint(
        ((RotorBore * 2 + (AirGap + MeshDensityAG / 2)) / 2) * np.sin(np.pi / 4),
        ((RotorBore * 2 + (AirGap + MeshDensityAG / 2)) / 2) * np.cos(np.pi / 4),
        0,
        MeshDensityAG,
        StatorStartPoint + 2,
    )
    gmsh.model.geo.addPoint(
        ((RotorBore * 2 + (AirGap + MeshDensityAG / 2)) / 2) * np.sin(-np.pi / 4),
        ((RotorBore * 2 + (AirGap + MeshDensityAG / 2)) / 2) * np.cos(-np.pi / 4),
        0,
        MeshDensityAG,
        StatorStartPoint + 3,
    )
    #stator teeth points
    n = StatorStartPoint + 4
    for a in range(spp * 3 + 1):
        if a == 0:
            gmsh.model.geo.addPoint(
                (StatorOD / 2 - BackIron)
                * np.sin(np.pi / 4 - (a * (np.pi / 2 / (spp * 3)))),
                (StatorOD / 2 - BackIron)
                * np.cos(np.pi / 4 - (a * (np.pi / 2 / (spp * 3)))),
                0,
                StatorDensity,
                n,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2 + Winding_division + 2)
                * np.sin(np.pi / 4 - (a * (np.pi / 2 / (spp * 3)))),
                (StatorID / 2 + Winding_division + 2)
                * np.cos(np.pi / 4 - (a * (np.pi / 2 / (spp * 3)))),
                0,
                StatorDensity,
                n + 1,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2 + 2) * np.sin(np.pi / 4 - (a * (np.pi / 2 / (spp * 3)))),
                (StatorID / 2 + 2) * np.cos(np.pi / 4 - (a * (np.pi / 2 / (spp * 3)))),
                0,
                MeshDensityAG * 2,
                n + 2,
            )
            gmsh.model.geo.addPoint(
                (StatorID) / 2 * np.sin(np.pi / 4 - (a * (np.pi / 2 / (spp * 3)))),
                (StatorID) / 2 * np.cos(np.pi / 4 - (a * (np.pi / 2 / (spp * 3)))),
                0,
                MeshDensityAG * 2,
                n + 3,
            )
        else:
            gmsh.model.geo.addPoint(
                (StatorOD / 2 - BackIron)
                * np.sin(np.pi / 4 - (a * (np.pi / 2 / (spp * 3))) + ToothAngleID),
                (StatorOD / 2 - BackIron)
                * np.cos(np.pi / 4 - (a * (np.pi / 2 / (spp * 3))) + ToothAngleID),
                0,
                StatorDensity,
                n,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2 + Winding_division + 2)
                * np.sin(np.pi / 4 - (a * (np.pi / 2 / (spp * 3))) + ToothAngleOD2),
                (StatorID / 2 + Winding_division + 2)
                * np.cos(np.pi / 4 - (a * (np.pi / 2 / (spp * 3))) + ToothAngleOD2),
                0,
                StatorDensity,
                n + 1,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2 + 2)
                * np.sin(np.pi / 4 - (a * (np.pi / 2 / (spp * 3))) + ToothAngleOD),
                (StatorID / 2 + 2)
                * np.cos(np.pi / 4 - (a * (np.pi / 2 / (spp * 3))) + ToothAngleOD),
                0,
                MeshDensityAG * 2,
                n + 2,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2 + 2)
                * np.sin(np.pi / 4 - ((a) * (np.pi / 2 / (spp * 3))) + teeth_angle_outer),
                (StatorID / 2 + 2)
                * np.cos(np.pi / 4 - ((a) * (np.pi / 2 / (spp * 3))) + teeth_angle_outer),
                0,
                MeshDensityAG * 2,
                n + 3,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2)
                * np.sin(np.pi / 4 - ((a) * (np.pi / 2 / (spp * 3))) + teeth_angle_outer),
                (StatorID / 2)
                * np.cos(np.pi / 4 - ((a) * (np.pi / 2 / (spp * 3))) + teeth_angle_outer),
                0,
                MeshDensityAG * 2,
                n + 4,
            )
        n += 5
        if a == spp * 3:
            gmsh.model.geo.addPoint(
                (StatorOD / 2 - BackIron)
                * np.sin(np.pi / 4 - a * (np.pi / 2 / (spp * 3))),
                (StatorOD / 2 - BackIron)
                * np.cos(np.pi / 4 - a * (np.pi / 2 / (spp * 3))),
                0,
                StatorDensity,
                n,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2 + Winding_division + 2)
                * np.sin(np.pi / 4 - a * (np.pi / 2 / (spp * 3))),
                (StatorID / 2 + Winding_division + 2)
                * np.cos(np.pi / 4 - a * (np.pi / 2 / (spp * 3))),
                0,
                StatorDensity,
                n + 1,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2 + 2) * np.sin(np.pi / 4 - a * (np.pi / 2 / (spp * 3))),
                (StatorID / 2 + 2) * np.cos(np.pi / 4 - a * (np.pi / 2 / (spp * 3))),
                0,
                MeshDensityAG * 2,
                n + 2,
            )
            gmsh.model.geo.addPoint(
                (StatorID) / 2 * np.sin(np.pi / 4 - a * (np.pi / 2 / (spp * 3))),
                (StatorID) / 2 * np.cos(np.pi / 4 - a * (np.pi / 2 / (spp * 3))),
                0,
                MeshDensityAG * 2,
                n + 3,
            )
        else:
            gmsh.model.geo.addPoint(
                (StatorOD / 2 - BackIron)
                * np.sin(np.pi / 4 - a * (np.pi / 2 / (spp * 3)) - ToothAngleID),
                (StatorOD / 2 - BackIron)
                * np.cos(np.pi / 4 - a * (np.pi / 2 / (spp * 3)) - ToothAngleID),
                0,
                StatorDensity,
                n,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2 + Winding_division + 2)
                * np.sin(np.pi / 4 - a * (np.pi / 2 / (spp * 3)) - ToothAngleOD2),
                (StatorID / 2 + Winding_division + 2)
                * np.cos(np.pi / 4 - a * (np.pi / 2 / (spp * 3)) - ToothAngleOD2),
                0,
                StatorDensity,
                n + 1,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2 + 2)
                * np.sin(np.pi / 4 - a * (np.pi / 2 / (spp * 3)) - ToothAngleOD),
                (StatorID / 2 + 2)
                * np.cos(np.pi / 4 - a * (np.pi / 2 / (spp * 3)) - ToothAngleOD),
                0,
                MeshDensityAG * 2,
                n + 2,
            )
            gmsh.model.geo.addPoint(
                (StatorID / 2 + 2)
                * np.sin(np.pi / 4 - a * (np.pi / 2 / (spp * 3)) - teeth_angle_inner),
                (StatorID / 2 + 2)
                * np.cos(np.pi / 4 - a * (np.pi / 2 / (spp * 3)) - teeth_angle_inner),
                0,
                MeshDensityAG * 2,
                n + 3,
            )
            gmsh.model.geo.addPoint(
                (StatorID) / 2
                * np.sin(np.pi / 4 - a * (np.pi / 2 / (spp * 3)) - teeth_angle_inner),
                (StatorID) / 2
                * np.cos(np.pi / 4 - a * (np.pi / 2 / (spp * 3)) - teeth_angle_inner),
                0,
                MeshDensityAG * 2,
                n + 4,
            )
        n += 5
    UpdatePosition = n
    gmsh.model.geo.synchronize()    #uploads points to model
    ######################################################################
    ######################### Model Lines ###############################
    ######################################################################
    gmsh.model.geo.addLine(2, 4, 1)
    gmsh.model.geo.addLine(4, 6, 3)
    gmsh.model.geo.addLine(3, 5, 2)
    gmsh.model.geo.addLine(5, 7, 4)
    gmsh.model.geo.addCircleArc(2, 1, 3, 7)
    gmsh.model.geo.addCircleArc(4, 1, 5, 8)
    gmsh.model.geo.addCircleArc(6, 1, 7, 9)
    gmsh.model.geo.addCircleArc(12, 1, 6, 12)  # 1
    gmsh.model.geo.addCircleArc(13, 1, 12, 13)  # 2
    gmsh.model.geo.addCircleArc(7, 1, 13, 18)  # 5

    degraded_curves = [7,8]
    n = 12 + n
    line = 21
    Loops = 1
    AirLoops = []
    RotorLoops = []
    run = []
    AirSurfaces = []

    # rotor points
    for a in range(len(Container2)):
        Container2[a].append([])
        if len(Container2[a][1]) > 1:
            looper = []
            for b in range(len(Container2[a][1]) - 1):
                Points = Features[Container2[a][1][b + 1]]
                run.append(Container2[a][1][b + 1])
                PTL = []  # points to link
                for c in range(len(Points)):
                    if (
                            c != 0
                            and round(Points[c][0], 4) == round(Points[0][0], 4)
                            and round(Points[c][1], 4) == round(Points[0][1], 4)
                    ):
                        pass
                    elif c < len(Points) - 2:
                        if round(Points[c][0], 4) == round(
                                Points[c + 1][0], 4
                        ) and round(Points[c][1], 4) == round(Points[c + 1][1], 4):
                            pass
                        else:
                            gmsh.model.geo.addPoint(
                                Points[c][0], Points[c][1], 0, MeshDensityR, n
                            )
                            PTL.append(n)
                            n += 1
                    else:
                        gmsh.model.geo.addPoint(
                            Points[c][0], Points[c][1], 0, MeshDensityR, n
                        )
                        PTL.append(n)
                        n += 1
                Loop = []
                for c in range(len(PTL) - 1):
                    gmsh.model.geo.addLine(PTL[c], PTL[c + 1], line)
                    Loop.append(line)
                    degraded_curves.append(line)
                    line += 1
                    if len(PTL) - 2 == c:
                        gmsh.model.geo.addLine(PTL[c + 1], PTL[0], line)
                        Loop.append(line)
                        degraded_curves.append(line)
                        line += 1
                gmsh.model.geo.addCurveLoop(Loop, Loops)
                looper.append(Loops)
                Surfaced = False
                for search in range(len(Container2)):
                    if Container2[search][1][0] == Container2[a][1][b + 1]:
                        storage = [Loops]
                        for c in Container2[search][2]:
                            storage.append(c)
                        Surfaced = True
                        gmsh.model.geo.addPlaneSurface(storage, 2000 + storage[0])
                if not Surfaced:
                    gmsh.model.geo.addPlaneSurface([Loops], 2000 + Loops)
                Container2[a][2].append(Loops)
                Loops += 1

    for a in range(len(Features)):
        Jump = False
        for b in run:
            if a == b:
                Jump = True
        if not Jump:
            Points = Features[a]
            PTL = []  # points to link
            for b in range(len(Points)):
                if (
                        round(Points[b][0], 6) == round(Points[0][0], 6)
                        and round(Points[b][1], 6) == round(Points[0][1], 6)
                        and b != 0
                ):
                    pass
                elif b < len(Points) - 1:
                    if round(Points[b][0], 6) == round(Points[b + 1][0], 6) and round(
                            Points[b][1], 6
                    ) == round(Points[b + 1][1], 6):
                        pass
                    else:
                        gmsh.model.geo.addPoint(
                            Points[b][0], Points[b][1], 0, MeshDensityR, n
                        )
                        PTL.append(n)
                        n += 1
                else:
                    gmsh.model.geo.addPoint(
                        Points[b][0], Points[b][1], 0, MeshDensityR, n
                    )
                    PTL.append(n)
                    n += 1
            Loop = []
            if len(PTL) > 3:
                for b in range(len(PTL) - 1):
                    gmsh.model.geo.addLine(PTL[b], PTL[b + 1], line)
                    Loop.append(line)
                    degraded_curves.append(line)
                    line += 1
                    if len(PTL) - 2 == b:
                        gmsh.model.geo.addLine(PTL[b + 1], PTL[0], line)
                        Loop.append(line)
                        degraded_curves.append(line)
                        line += 1
                gmsh.model.geo.addCurveLoop(Loop, Loops)
                Surface = [Loops]
                for b in range(len(Container2)):
                    if a == Container2[b][1][0]:
                        for c in Container2[b][2]:
                            Surface.append(c)
                gmsh.model.geo.addPlaneSurface(Surface, 2000 + Loops)
                AirLoops.append(Loops)
                AirSurfaces.append(2000 + Loops)
                Loops += 1

    ################################################################
    ########################### Stator Lines #######################
    ################################################################
    StatorSurface = []
    Conductor = []
    AirGapSurface = []
    IP = StatorStartPoint + 4
    StatorStartLine = n
    gmsh.model.geo.addLine(StatorStartPoint, IP, n)
    n += 1
    gmsh.model.geo.addLine(IP, IP + 1, n)
    n += 1
    gmsh.model.geo.addLine(IP + 1, IP + 2, n)
    n += 1
    gmsh.model.geo.addLine(IP + 2, IP + 3, n)
    n += 1
    StatorSurface.append(StatorStartLine)
    StatorSurface.append(StatorStartLine + 1)
    StatorSurface.append(StatorStartLine + 2)
    StatorSurface.append(StatorStartLine + 3)
    Stator_boundaries_Right = [StatorStartLine + 3, StatorStartLine + 2, StatorStartLine + 1, StatorStartLine]
    IP = StatorStartPoint + 4 + 5
    for a in range(spp * 3):
        if a == 0:
            gmsh.model.geo.addCircleArc(IP - 2 + (8 * a), 1, IP + 4 + (8 * a), n)  # stator surface
            AirGapSurface.append(n)
            StatorSurface.append(n)
            n += 1

        if a == spp * 3 - 1:
            gmsh.model.geo.addLine(IP + 4 + (10 * a), IP + 3 + (10 * a), n)  # tooth cap
            StatorSurface.append(n)
            AirGapSurface.append(n)
            gmsh.model.geo.addLine(IP + 3 + (10 * a), IP + 2 + (10 * a), n + 1)  # tooth cap
            StatorSurface.append(n + 1)
            gmsh.model.geo.addLine(IP + 2 + (10 * a), IP + 1 + (10 * a), n + 2)  # C1
            StatorSurface.append(n + 2)
            gmsh.model.geo.addLine(IP + 1 + (10 * a), IP + (10 * a), n + 3)  # C2
            StatorSurface.append(n + 3)
            gmsh.model.geo.addLine(IP + (10 * a), IP + 5 + (10 * a), n + 4)  # C2 Top
            StatorSurface.append(n + 4)

            gmsh.model.geo.addLine(IP + 6 + (10 * a), IP + 1 + (10 * a), n + 10)  # C1 Top
            gmsh.model.geo.addLine(IP + 3 + (10 * a), IP + 8 + (10 * a), n + 11)  # C1 Bottom
            AirGapSurface.append(n + 11)

            gmsh.model.geo.addLine(IP + 5 + (10 * a), IP + 6 + (10 * a), n + 5)  # C2
            StatorSurface.append(n + 5)
            gmsh.model.geo.addLine(IP + 6 + (10 * a), IP + 7 + (10 * a), n + 6)  # C1
            StatorSurface.append(n + 6)
            gmsh.model.geo.addLine(IP + 7 + (10 * a), IP + 8 + (10 * a), n + 7)  # tooth cap
            StatorSurface.append(n + 7)
            gmsh.model.geo.addLine(IP + 8 + (10 * a), IP + 9 + (10 * a), n + 8)  # tooth cap
            StatorSurface.append(n + 8)
            AirGapSurface.append(n + 8)
            gmsh.model.geo.addCircleArc(IP + 9 + (10 * a), 1, IP + 9 + 4 + (10 * a), n + 9)  # stator surface
            StatorSurface.append(n + 9)
            AirGapSurface.append(n + 9)

            Conductor.append([n + 3, n + 4, n + 5, n + 10])
            Conductor.append([n + 2, -(n + 10), n + 6, n + 7, -(n + 11), n + 1])
            n += 12
        else:
            gmsh.model.geo.addLine(IP + 4 + (10 * a), IP + 3 + (10 * a), n)  # tooth cap
            StatorSurface.append(n)
            AirGapSurface.append(n)
            gmsh.model.geo.addLine(IP + 3 + (10 * a), IP + 2 + (10 * a), n + 1)  # tooth cap
            StatorSurface.append(n + 1)
            gmsh.model.geo.addLine(IP + 2 + (10 * a), IP + 1 + (10 * a), n + 2)  # C1
            StatorSurface.append(n + 2)
            gmsh.model.geo.addLine(IP + 1 + (10 * a), IP + (10 * a), n + 3)  # C2
            StatorSurface.append(n + 3)
            gmsh.model.geo.addLine(IP + (10 * a), IP + 5 + (10 * a), n + 4)  # C1 Top
            StatorSurface.append(n + 4)

            gmsh.model.geo.addLine(IP + 6 + (10 * a), IP + 1 + (10 * a), n + 10)  # C2Top
            gmsh.model.geo.addLine(IP + 3 + (10 * a), IP + 8 + (10 * a), n + 11)  # C1 Bottom
            AirGapSurface.append(n + 11)

            gmsh.model.geo.addLine(IP + 5 + (10 * a), IP + 6 + (10 * a), n + 5)  # C2
            StatorSurface.append(n + 5)
            gmsh.model.geo.addLine(IP + 6 + (10 * a), IP + 7 + (10 * a), n + 6)  # C1
            StatorSurface.append(n + 6)
            gmsh.model.geo.addLine(IP + 7 + (10 * a), IP + 8 + (10 * a), n + 7)  # tooth cap
            StatorSurface.append(n + 7)
            gmsh.model.geo.addLine(IP + 8 + (10 * a), IP + 9 + (10 * a), n + 8)  # tooth cap
            StatorSurface.append(n + 8)
            AirGapSurface.append(n + 8)
            gmsh.model.geo.addCircleArc(IP + 9 + (10 * a), 1, IP + 9 + 5 + (10 * a), n + 9)  # stator surface
            StatorSurface.append(n + 9)
            AirGapSurface.append(n + 9)

            Conductor.append([n + 3, n + 4, n + 5, n + 10])
            Conductor.append([n + 2, -(n + 10), n + 6, n + 7, -(n + 11), n + 1])
            n += 12

    gmsh.model.geo.addLine(UpdatePosition - 5, StatorStartPoint + 1, n)
    gmsh.model.geo.addLine(UpdatePosition - 4, UpdatePosition - 5, n + 1)
    gmsh.model.geo.addLine(UpdatePosition - 3, UpdatePosition - 4, n + 2)
    gmsh.model.geo.addLine(UpdatePosition - 2, UpdatePosition - 3, n + 3)
    StatorSurface.append(n + 3)
    StatorSurface.append(n + 2)
    StatorSurface.append(n + 1)
    StatorSurface.append(n)
    Stator_boundaries_Left = [n + 3, n + 2, n + 1, n]
    gmsh.model.geo.addCircleArc(StatorStartPoint + 1, 1, StatorStartPoint, n + 4)
    StatorSurface.append(n + 4)
    StatorSurfaceCurve = n + 4

    # upper air boundary
    gmsh.model.geo.addLine(IP - 2, StatorStartPoint + 2, n + 5)
    gmsh.model.geo.addCircleArc(StatorStartPoint + 2, 1, StatorStartPoint + 3, n + 6)
    gmsh.model.geo.addLine(StatorStartPoint + 3, UpdatePosition - 2, n + 7)
    AirGapSurface.append(-(n + 7))
    AirGapSurface.append(-(n + 6))
    AirGapSurface.append(-(n + 5))
    #gmsh.model.geo.synchronize()
    #gmsh.fltk.run()

    ###########################################################################
    ################################# Add Curve Loops #########################
    ###########################################################################
    Rotor = Loops
    gmsh.model.geo.addCurveLoop([1, 8, -2, -7], Loops)
    LowwerAirgap = Loops + 1
    gmsh.model.geo.addCurveLoop([3, 9, -4, -8], Loops + 1)
    gmsh.model.geo.addCurveLoop([5, 11, -6, -10], Loops + 2)

    RotorSurface = [Rotor]
    for a in AirLoops:
        skip = True
        for b in RotorLoops:
            if a == b:
                skip = False
        if skip:
            RotorSurface.append(a)
    StatorSurfaceStart = n
    n += 1
    gmsh.model.geo.addCurveLoop(AirGapSurface, n)
    n += 1
    gmsh.model.geo.addCurveLoop(StatorSurface, n)
    n += 1
    for a in range(len(Conductor)):
        gmsh.model.geo.addCurveLoop(Conductor[a], n)
        n += 1

    ##############################################################
    ####################### Add Surfaces #########################
    ##############################################################
    gmsh.model.geo.addPlaneSurface([StatorSurfaceStart + 1], 20000)
    gmsh.model.geo.addPlaneSurface([StatorSurfaceStart + 2], 20001)
    n = StatorSurfaceStart + 3
    m = 2
    for a in range(len(Conductor)):
        gmsh.model.geo.addPlaneSurface([n], 20000 + m)
        n += 1
        m += 1
    gmsh.model.geo.addPlaneSurface(RotorSurface, 1000)
    gmsh.model.geo.addPlaneSurface([LowwerAirgap], 1001)
    RotorSteel = [1000]
    for a in Container2:
        for b in a[2]:
            if a[0] / 2 == np.floor(a[0] / 2):
                RotorSteel.append(b + 2000)
            else:
                AirSurfaces.append(b + 2000)

    for a in StatorSurface:
        degraded_curves.append(a)

    ####################################################
    ################## Add Groups ######################
    ####################################################
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(2, AirSurfaces, 10000, name="air")          # rotor air
    gmsh.model.addPhysicalGroup(2, RotorSteel, 10004, name="steel")         # rotor
    gmsh.model.addPhysicalGroup(2, [1001], 10002, name="air_Bottom")        # rotor airgap
    gmsh.model.addPhysicalGroup(2, [20001], tag=10005, name="Stator_Iron")  # Stator
    gmsh.model.addPhysicalGroup(2, [20000], tag=10006, name="Air_Gap")      # Stator airGap
    gmsh.model.addPhysicalGroup(2, [20002], tag=10008, name="A1")           # Conductor A
    gmsh.model.addPhysicalGroup(2, [20003], tag=10009, name="A2")           # Conductor A
    gmsh.model.addPhysicalGroup(2, [20004], tag=10010, name="A3")           # Conductor A
    gmsh.model.addPhysicalGroup(2, [20005], tag=10011, name="A4")           # Conductor A
    gmsh.model.addPhysicalGroup(2, [20006], tag=10012, name="A5")           # Conductor A
    gmsh.model.addPhysicalGroup(2, [20007], tag=10013, name="A6")           # Conductor A
    gmsh.model.addPhysicalGroup(2, [20008], tag=10016, name="B1")           # Conductor B
    gmsh.model.addPhysicalGroup(2, [20009], tag=10017, name="B2")           # Conductor B
    gmsh.model.addPhysicalGroup(2, [20010], tag=10018, name="B3")           # Conductor B
    gmsh.model.addPhysicalGroup(2, [20011], tag=10019, name="B4")           # Conductor B
    gmsh.model.addPhysicalGroup(2, [20012], tag=10020, name="B5")           # Conductor B
    gmsh.model.addPhysicalGroup(2, [20013], tag=10021, name="B6")           # Conductor B
    gmsh.model.addPhysicalGroup(2, [20014], tag=10024, name="C1")           # Conductor C
    gmsh.model.addPhysicalGroup(2, [20015], tag=10025, name="C2")           # Conductor C
    gmsh.model.addPhysicalGroup(2, [20016], tag=10026, name="C3")           # Conductor C
    gmsh.model.addPhysicalGroup(2, [20017], tag=10027, name="C4")           # Conductor C
    gmsh.model.addPhysicalGroup(2, [20018], tag=10028, name="C5")           # Conductor C
    gmsh.model.addPhysicalGroup(2, [20019], tag=10029, name="C6")           # Conductor C

    # surfaces
    gmsh.model.addPhysicalGroup(1, [9], 20000, name="RotorSurface_0")                           # rotor moving band
    gmsh.model.addPhysicalGroup(1, [12], 20003, name="RotorSurface_1")                          # rotor moving band
    gmsh.model.addPhysicalGroup(1, [13], 20002, name="RotorSurface_2")                          # rotor moving band
    gmsh.model.addPhysicalGroup(1, [18], 20001, name="RotorSurface_3")                          # rotor moving band
    gmsh.model.addPhysicalGroup(1, [7], 20016, name="Dirchlet")                                 # shaft boundary
    gmsh.model.addPhysicalGroup(1, [StatorSurfaceCurve + 2], tag=20019, name="Stator_Boundry")  # stator airgap
    gmsh.model.addPhysicalGroup(1, [StatorSurfaceCurve], tag=20026, name="Dirchlet_1")          # stator boundary
    #periodic boundaries
    gmsh.model.addPhysicalGroup(1, [1], 20008, name="PeriodicBC_R1")                            #rotor right
    gmsh.model.addPhysicalGroup(1, [3], 20009, name="PeriodicBC_R2")                            #rotor airgap right
    gmsh.model.addPhysicalGroup(1, [2], 20012, name="PeriodicBC_L1")                            #rotor left
    gmsh.model.addPhysicalGroup(1, [4], 20013, name="PeriodicBC_L2")                            #rotor airgap left
    gmsh.model.addPhysicalGroup(1, degraded_curves, 20017, name="Damage")                            #rotor airgap left
    gmsh.model.addPhysicalGroup(1, [StatorSurfaceCurve + 1], tag=20020, name="Periodic_Boundry_Air_R")
    gmsh.model.addPhysicalGroup(1, [Stator_boundaries_Right[0]], tag=20021, name="Periodic_Boundry_Stator_Lower_R")
    gmsh.model.addPhysicalGroup(1, [Stator_boundaries_Right[1]], tag=20022, name="Periodic_Boundry_Stator_Upper_R")
    gmsh.model.addPhysicalGroup(1, [StatorSurfaceCurve + 3], tag=20023, name="Periodic_Boundry_Air_L")
    gmsh.model.addPhysicalGroup(1, [Stator_boundaries_Left[0]], tag=20024, name="Periodic_Boundry_Stator_Lower_L")
    gmsh.model.addPhysicalGroup(1, [Stator_boundaries_Left[1]], tag=20025, name="Periodic_Boundry_Stator_Upper_L")
    gmsh.model.addPhysicalGroup(1, [Stator_boundaries_Right[2]], tag=20027, name="Periodic_Boundry_Stator_Upper_R2")
    gmsh.model.addPhysicalGroup(1, [Stator_boundaries_Right[3]], tag=20028, name="Periodic_Boundry_Stator_Upper_R3")
    gmsh.model.addPhysicalGroup(1, [Stator_boundaries_Left[2]], tag=20029, name="Periodic_Boundry_Stator_Upper_L2")
    gmsh.model.addPhysicalGroup(1, [Stator_boundaries_Left[3]], tag=20030, name="Periodic_Boundry_Stator_Upper_L3")

    # Points
    gmsh.model.addPhysicalGroup(0, [1], 30002, name="Center of Machine")

    ######################################################
    ##################### Runing Model ###################
    ######################################################
    try:
        gmsh.model.mesh.generate(2)
        gmsh.plugin.setNumber('Distance', 'PhysicalLine', 20017)  # value applied to airgap
        gmsh.plugin.setNumber('Distance', 'DistanceType', 0)
        gmsh.plugin.run('Distance')
        gmsh.view.write(gmsh.view.get_tags()[0], directory + "\\Distance.pos")
        Function.minimise_msh(directory, "Distance", mag, ".POS")
    except:
        gmsh.write("Models\\LSModel_Error.msh")
        error("Meshing Error")

    gmsh.write(directory + "\\LSModel.msh")
    Function.minimise_msh(directory, "LSModel", mag, ".msh")
    FileLoadFile = open(
        str(directory) + "Angle.PRO", "w"
    )
    FileLoadFile.write("DefineConstant [Angle = " + str(Angle1) + "];\n")
    FileLoadFile.write("DefineConstant [AngleSTR = '" + str(Angle1) + "'];\n")
    FileLoadFile.write("current_multiplier = " + str(current) + ";\n")
    FileLoadFile.close()

    try:
        gmsh.onelab.run(
            "GetDp",
            '"'+getDP+'" '+
            '"' + str(directory) + 'LSModel_mini.PRO" -solve -pos',
        )
    except:
        gmsh.model.mesh.generate(2)
        gmsh.fltk.run()
        error("Error in the model please check and try again")

    print("\rRotating Rotor", end="")
    FileLoadFile = open(
        str(directory) + "Angle.PRO", "w"
    )
    FileLoadFile.write("DefineConstant [Angle = " + str(Angle2) + "];\n")
    FileLoadFile.write("DefineConstant [AngleSTR = '" + str(Angle2) + "'];\n")
    FileLoadFile.write("current_multiplier = " + str(current) + ";\n")
    FileLoadFile.close()
    try:
        gmsh.onelab.run(
            "GetDp",
            '"'+getDP+'" '+
            '"' + str(directory) + 'LSModel_mini.PRO" -solve -pos',
        )
    except:
        error("Error in the model please check and try again")
    gmsh.finalize()
    print("\rFinnished modeling", end="")
