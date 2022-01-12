
# Run with: blender --background --python convert_directory.py 1> /dev/null
# Note: " 1> " supresses all output from the file to enable appropriate presentation of the tqdm progress bar (-> no printing possible)

import os
import sys
import bpy

sys.path.append('/home/tino/miniconda3/lib/python3.9/site-packages/')

try:
    from tqdm import tqdm
except:
    print('WARNING: tqdm is not available; be patient this can take several minutes', file=sys.stderr)
    tqdm = lambda x: x

path = './Meshes/apfelsine/apfelsine-meshes/'


# Create subfolder for export
export_path = path + 'blender_export/'
try:
    os.mkdir(export_path)
except FileExistsError:
    pass


print("Converting Data", file=sys.stderr)

for file in tqdm([f for f in os.listdir(path) if f.endswith('.obj')]):
        obj_in = path + file
        obj_out = export_path + file

        bpy.ops.wm.read_factory_settings(use_empty=True)

        bpy.ops.import_scene.obj(filepath=obj_in, axis_forward='-Z', axis_up='Y')
        bpy.ops.export_scene.obj(filepath=obj_out, axis_forward='-Z', axis_up='Y')

