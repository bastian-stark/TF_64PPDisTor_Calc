#C:\Users\basti\AppData\Local\Programs\Python\Python39\Structural_Files\structural_analysis_calculations\Human\TFs\NFYA
#res cut-off: 3.50

import os
import math
import numpy as np

print('Current working directory: ', os.getcwd())
directory = os.getcwd()
directory_reply = input('Change directory? Y/N: ')
directory_reply = directory_reply.upper()
if directory_reply == "Y":
    directory = input('Change directory to: ')
    os.chdir(directory)
    print('New working directory: ', os.getcwd())
else:
    pass

format_reply = input("Reformat PDB files? Y/N: ")
format_reply = format_reply.upper()
if format_reply == "Y":
    print("Re-formatting PDB files into TSV files with un-merged columns...")
    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
            split_basename = os.path.splitext(filename)
            basename = split_basename[0]
            print("Reformatting File: ", filename)
            PDBres3p50FilePath = str(basename) + ".tsv"
            with open(filename, 'r') as pdbFile:
                with open(PDBres3p50FilePath, 'w') as pdbres3p50File:
                    for line in pdbFile:
                        if line[0:6].strip() != "ATOM": continue
                        name = line[12:16].strip()
                        residueName = line[17:20].strip()
                        chainID = line[21].strip()
                        resSeqNum = line[22:26].strip()
                        xPos = line[30:38].strip()
                        yPos = line[38:46].strip()
                        zPos = line[46:54].strip()
                        tempFactor = line[60:66].strip()
                        relevantData = (name, residueName, chainID, resSeqNum, xPos, yPos, zPos, tempFactor)
                        for data in relevantData: assert data, "Data point missing: " + '\n\t' + str(
                            relevantData) + '\n\t' + line
                        pdbres3p50File.write('\t'.join(relevantData) + '\n')
else:
    pass

res_input = input("Enter resolution cut-off (n.nn): ")
PP_type_response = input("Include 6-4PP analysis? (Y/N)")
PP_type_response = PP_type_response.upper()
if PP_type_response == "Y":
    PP_type2 = 1
else: PP_type2 = 0

print("Retrieving atom coordinates...")
################################################################CPD#####################################################
if PP_type2 == 1:
    PP64_atom_list = []
    for filename in os.listdir(directory):
        if filename.endswith(".tsv"):
            print("Retreiving From: ", filename)
            offset_file = open("motif_offsets.txt")
            for line in offset_file:
                offset_list = line.split()
                structure_id = offset_list[0]
                chain = offset_list[1]
                central_c = offset_list[2]
                resolution = offset_list[3]
                if filename.startswith(structure_id) and resolution <= res_input:
                    with open(filename) as PDBfile:
                        for line in PDBfile:
                            pdb_line = line.split()
                            if pdb_line[0] == "C4" or pdb_line[0] == "C5" or pdb_line[0] == "C6" or pdb_line[0] == "O4" or pdb_line[0] == "N4":
                                if pdb_line[1] == "DT" or pdb_line[1] == "DC":
                                    if pdb_line[2] == chain:
                                        i = -11
                                        while i <= 11:
                                            position = str(int(central_c) + i)
                                            position_relative_to_midpoint = i
                                            if pdb_line[3] == position and float(pdb_line[7]) > 0:
                                                pdb_line[3] = position_relative_to_midpoint
                                                position_rtm = pdb_line[3]
                                                PP64_atom_list.append(pdb_line)
                                            i = i + 1
            offset_file.close()
    print("Pairing atoms...")
    atom_pair_coord_list = []
    previous_line = ['A1', 'DN', 'Z', '999', '999', '999', '999', '999']
    for line in PP64_atom_list:
        current_line = line
        if int(current_line[3]) == int(previous_line[3]):
            if previous_line[0] == "C5" and current_line[0] == "C6":
                C5C6 = previous_line[4:7], current_line[4:7], previous_line[3], 'C5C6', current_line[1]
                atom_pair_coord_list.append(C5C6)
            elif previous_line[0] == "C4" and current_line[0] == "N4":
                C4N4 = previous_line[4:7], current_line[4:7], previous_line[3], 'C4N4', current_line[1]
                atom_pair_coord_list.append(C4N4)
            elif previous_line[0] == "C4" and current_line[0] == "O4":
                C4O4 = previous_line[4:7], current_line[4:7], previous_line[3], 'C4O4', current_line[1]
                atom_pair_coord_list.append(C4O4)
        previous_line = current_line
    print("Pairing pyrimidines...")
    dipy_list = []
    previous_line = (['999', '999', '999'], ['999', '999', '999'], 999, 'AZAZ', 'DN')
    for line in atom_pair_coord_list:
        current_line = line
        if abs(int(previous_line[2]) - int(current_line[2])) == 1 and previous_line[3] == "C5C6" and (current_line[3] == "C4O4" or current_line[3] == "C4N4"):
            previous_x1_average = (float(previous_line[0][0]) + float(previous_line[1][0])) / 2
            previous_y1_average = (float(previous_line[0][1]) + float(previous_line[1][1])) / 2
            previous_z1_average = (float(previous_line[0][2]) + float(previous_line[1][2])) / 2
            previous_x1y1z1_average = [previous_x1_average, previous_y1_average, previous_z1_average]
            current_x1_average = (float(current_line[0][0]) + float(current_line[1][0])) / 2
            current_y1_average = (float(current_line[0][1]) + float(current_line[1][1])) / 2
            current_z1_average = (float(current_line[0][2]) + float(current_line[1][2])) / 2
            current_x1y1z1_average = [current_x1_average, current_y1_average, current_z1_average]
            position_relative_to_midpoint = (int(previous_line[2]) + int(current_line[2])) / 2
            dipy_info = previous_x1y1z1_average, current_x1y1z1_average, previous_line[0:2], current_line[
                                                                                             0:2], position_relative_to_midpoint
            dipy_list.append(dipy_info)
        previous_line = current_line
    print("Calculating distance and torsion angle...")
    dipy_distor_list = []
    dipy_distor = []
    lowdis_lowtor = []
    lowdis_hightor = []
    highdis_lowtor = []
    highdis_hightor = []
    for line in dipy_list:
        x1 = line[2][0][0]
        y1 = line[2][0][1]
        z1 = line[2][0][2]
        x2 = line[2][1][0]
        y2 = line[2][1][1]
        z2 = line[2][1][2]
        x3 = line[3][0][0]
        y3 = line[3][0][1]
        z3 = line[3][0][2]
        x4 = line[3][1][0]
        y4 = line[3][1][1]
        z4 = line[3][1][2]
        position_rtm = line[4]
        p = np.array([
            [float(x1), float(y1), float(z1)],
            [float(x2), float(y2), float(z2)],
            [float(x3), float(y3), float(z3)],
            [float(x4), float(y4), float(z4)]
        ])
        b = p[:-1] - p[1:]
        b[0] *= -1
        v = np.array([v - (v.dot(b[1]) / b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])
        v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1, 1)
        b1 = b[1] / np.linalg.norm(b[1])
        x = np.dot(v[0], v[1])
        m = np.cross(v[0], b1)
        y = np.dot(m, v[1])
        angle = np.degrees(np.arctan2(y, x))
        x1_average, y1_average, z1_average, x2_average, y2_average, z2_average = line[0][0], line[0][1], line[0][2], \
                                                                                 line[1][0], line[1][1], line[1][2]
        distance = math.sqrt(
            math.pow(float(x2_average) - float(x1_average), 2) + math.pow(float(y2_average) - float(y1_average),
                                                                          2) + math.pow(
                float(z2_average) - float(z1_average), 2))
        res3p50_line = position_rtm, distance, angle
        dipy_distor_list.append(res3p50_line)
        CPD_list_human = [(0.5, 1.203038), (1.5, 1.1469525), (2.5, 1.052703), (3.5, 0.96022575), (4.5, 0.91849755),
                          (5.5, 0.9589802), (6.5, 1.0136394), (7.5, 1.0588725), (8.5, 1.1400015), (9.5, 1.199938),
                          (10.5, 1.1835005), (11.5, 1.1289195), (12.5, 1.05815), (13.5, 0.98928), (14.5, 0.9479745),
                          (15.5, 0.9313637), (16.5, 0.9487725), (17.5, 1.0307363), (18.5, 1.143613), (19.5, 1.2095185),
                          (20.5, 1.2285635), (21.5, 1.2025635), (22.5, 1.1201735), (23.5, 1.0317885),
                          (24.5, 0.97762955),
                          (25.5, 0.95501775), (26.5, 0.96511645), (27.5, 1.02646475), (28.5, 1.121538),
                          (29.5, 1.198515),
                          (30.5, 1.211097), (31.5, 1.1576815), (32.5, 1.0812075), (33.5, 0.99849965), (34.5, 0.9371397),
                          (35.5, 0.92547065), (36.5, 0.9464726), (37.5, 0.9948185), (38.5, 1.0747465), (39.5, 1.16335),
                          (40.5, 1.204351), (41.5, 1.17833), (42.5, 1.1049375), (43.5, 1.0137727), (44.5, 0.9460494),
                          (45.5, 0.9153289), (46.5, 0.92186725), (47.5, 0.96024335), (48.5, 1.0291023),
                          (49.5, 1.1150975),
                          (50.5, 1.1659725), (51.5, 1.15324), (52.5, 1.0978995), (53.5, 1.02521005), (54.5, 0.9599301),
                          (55.5, 0.92540675), (56.5, 0.92970075), (57.5, 0.963809), (58.5, 1.01587245), (59.5, 1.07304),
                          (60.5, 1.1102475), (61.5, 1.11601), (62.5, 1.08703), (63.5, 1.0225974), (64.5, 0.96157145),
                          (65.5, 0.9360816), (66.5, 0.94626545), (67.5, 0.9855354), (68.5, 1.0307925), (69.5, 1.063308),
                          (70.5, 1.0834105), (71.5, 1.0788175), (72.5, 1.0490035), (73.5, 1.03141)]
        for line in CPD_list_human:
            info = position_rtm, distance, angle, line[1]
        dipy_distor.append(info)
        if distance < 4.17 and angle < 30.36:
            lowdis_lowtor.append(info)
        elif distance < 4.17 and angle > 41.72:
            lowdis_hightor.append(info)
        elif distance > 4.71 and angle < 30.36:
            highdis_lowtor.append(info)
        elif distance > 4.71 and angle > 41.72:
            highdis_hightor.append(info)
    dipy_distor_list.sort()
    dipy_distor.sort()
    lowdis_lowtor.sort()
    lowdis_hightor.sort()
    highdis_lowtor.sort()
    highdis_hightor.sort()
    with open("NFY_64PP_distance_angle_all_res3p50.txt", "w+") as f:
        f.write("Distance" + '\n')
        i = -10.5
        while i <= 10.5:
            f.write(str(i) + '\t')
            for line in dipy_distor:
                if line[0] == i:
                    f.write(str(line[1]) + '\t')
            f.write('\n')
            i = i + 1
        f.write("Torsion Angle" + '\n')
        i = -10.5
        while i <= 10.5:
            f.write(str(i) + '\t')
            for line in dipy_distor:
                if line[0] == i:
                    f.write(str(line[2]) + '\t')
            f.write('\n')
            i = i + 1
    f.close()
    with open("NFY_64PP_distance_angle_all_res3p50_BaWP.txt", "w+") as f:
        for line in dipy_distor:
            f.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str(line[2]) + '\n')
    f.close()
    with open("NFY_64PP_distance_angle_all_res3p50_cats.txt", "w+") as f:
        f.write("Low Distance Low Angle" + '\n')
        for line in lowdis_lowtor:
            f.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str(line[2]) + '\t' + str(line[3]) + '\n')
        f.write("Low Distance High Angle" + '\n')
        for line in lowdis_hightor:
            f.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str(line[2]) + '\t' + str(line[3]) + '\n')
        f.write("High Distance Low Angle" + '\n')
        for line in highdis_lowtor:
            f.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str(line[2]) + '\t' + str(line[3]) + '\n')
        f.write("High Distance High Angle" + '\n')
        for line in highdis_hightor:
            f.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str(line[2]) + '\t' + str(line[3]) + '\n')
    f.close()