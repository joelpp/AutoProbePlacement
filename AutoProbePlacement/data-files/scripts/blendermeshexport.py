import bpy
import os

# get the current path and make a new folder for the exported meshes
path = bpy.path.abspath('//objs\\')
os.makedirs(path)

for object in bpy.context.selected_objects:

    # deselect all meshes
    bpy.ops.object.select_all(action='DESELECT')

    # select the object
    object.select = True

    # export object with its name as file name
    bpy.ops.export_scene.obj(filepath=str((path + object.name + '.obj')), use_selection=True)