import bpy
import csv
import numpy as np
import time
import math

#### Importing CSV files and saving to response vectors, one for each degree of freedom
fpath = 'C:/7GDL_RespostaAnim2.csv'
csvFile = csv.reader(open(fpath))
data = [row for row in csvFile][0:]
data_7gdl = []
response = []

for i in range(len(data)):
    data_list = data[i]
    data_str = data_list[0]
    data_float = [float(x) for x in data_str.split(' ')]
    if i == 0:
        gdls = len(data_float)
    data_7gdl.append(data_float)
    
for i in range(gdls):
    response_col = []
    for j in range(len(data)):
        data_point = data_7gdl[j][i]
        response_col.append(data_point)
    response.append(response_col)

resp_1 = response[0]
resp_2 = response[1]
resp_3 = response[2]
resp_4 = response[3]
resp_5 = response[4]
resp_6 = response[5] 
resp_7 = response[6]


#### Location of the bodies in space
L1 = (3/5)*1.7
L2 = 1.7-L1
L3 = 1.4 
#### Creating mesh bodies, onde cylinder for each whell and a cube for the chassi
op1 = bpy.ops.mesh.primitive_cylinder_add(radius= 0.25, depth=0.15, enter_editmode=False, align='CURSOR', location=(L1, L3/2, 0), rotation=(math.pi/2, 0, 0), scale=(1, 1, 1))
op2 = bpy.ops.mesh.primitive_cylinder_add(radius= 0.25, depth=0.15, enter_editmode=False, align='CURSOR', location=(L1, -L3/2, 0), rotation=(math.pi/2, 0, 0), scale=(1, 1, 1))
op3 = bpy.ops.mesh.primitive_cylinder_add(radius= 0.25, depth=0.15, enter_editmode=False, align='CURSOR', location=(-L2, -L3/2, 0), rotation=(math.pi/2, 0, 0), scale=(1, 1, 1))
op4 = bpy.ops.mesh.primitive_cylinder_add(radius= 0.25, depth=0.15, enter_editmode=False, align='CURSOR', location=(-L2, L3/2, 0), rotation=(math.pi/2, 0, 0), scale=(1, 1, 1))
op5 = bpy.ops.mesh.primitive_cube_add(size=0.5, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))


#### storing the objects into variables for loop call
cylinder1 = bpy.data.collections[0].objects['Cylinder']
cylinder2 = bpy.data.collections[0].objects['Cylinder.001']
cylinder3 = bpy.data.collections[0].objects['Cylinder.002']
cylinder4 = bpy.data.collections[0].objects['Cylinder.003']
cube5 = bpy.data.collections[0].objects['Cube']
frame_num = 0
#### assembly of animation
for i in range(len(data)):
    
    bpy.context.scene.frame_set(frame_num)
    cylinder1.location = (L1, L3/2,resp_1[i]) #Changing Cylinder1's location according with it's free vibration response vector
    cylinder1.keyframe_insert(data_path = "location", frame = frame_num) #Keyframe tracker for cylinder1's translation
    cylinder2.location = (L1, -L3/2,resp_2[i]) #Changing Cylinder2's location according with it's free vibration response vector
    cylinder2.keyframe_insert(data_path = "location", frame = frame_num) #Keyframe tracker for cylinder2's translation
    cylinder3.location = (-L2, -L3/2,resp_3[i]) #Changing Cylinder3's location according with it's free vibration response vector
    cylinder3.keyframe_insert(data_path = "location", frame = frame_num) #Keyframe tracker for cylinder3's translation
    cylinder4.location = (-L2, L3/2,resp_4[i]) #Changing Cylinder4's location according with it's free vibration response vector
    cylinder4.keyframe_insert(data_path = "location", frame = frame_num) #Keyframe tracker for cylinder4's translation
    cube5.location = (0,0,L1/3 + resp_5[i]) #Changing Cubes's location according with it's free vibration response vector
    #Changing Cylinder1's orientation according with it's free vibration response vector
    bpy.ops.transform.rotate(value=resp_6[i], orient_axis='Y', orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False, release_confirm=True)
    bpy.ops.transform.rotate(value=resp_7[i], orient_axis='X', orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False, release_confirm=True)
    #Keyframe tracker for cubes's rotation
    cube5.keyframe_insert(data_path = "rotation_euler", frame = frame_num)
    #Keyframe tracker for cubes's translation
    cube5.keyframe_insert(data_path = "location", frame = frame_num)
    frame_num += 5
    

