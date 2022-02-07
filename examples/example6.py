import ctypes
import os
import sys
import shutil


# load vocaltractlab binary
PREFIX = 'lib'
SUFFIX = ''
if sys.platform.startswith('linux'):
    SUFFIX = '.so'
elif sys.platform.startswith('win32'):
    PREFIX = ''
    SUFFIX = '.dll'
elif sys.platform.startswith('darwin'):
    SUFFIX = '.dylib'

VTL = ctypes.cdll.LoadLibrary(f'../lib/Release/{PREFIX}VocalTractLabApi{SUFFIX}')
del PREFIX, SUFFIX


# get version / compile date
version = ctypes.c_char_p(b' ' * 64)
VTL.vtlGetVersion(version)
print('Compile date of the library: "%s"' % version.value.decode())


# initialize vtl
speaker_file_name = ctypes.c_char_p('../resources/JD3.speaker'.encode())

failure = VTL.vtlInitialize(speaker_file_name)
if failure != 0:
    raise ValueError('Error in vtlInitialize! Errorcode: %i' % failure)


# NOTE: You need the apfelsine.ges file, which will be create by running
# example5.py

print("Part 1: Synthesize gestural score as reference")
speaker_file_name = ctypes.c_char_p('../resources/JD3.speaker'.encode())
segment_file_name = ctypes.c_char_p(b'apfelsine.seg')
gesture_file_name = ctypes.c_char_p('apfelsine.ges'.encode())
wav_file_name = ctypes.c_char_p('apfelsine.wav'.encode())
feedback_file_name = ctypes.c_char_p('apfelsine.txt'.encode())

VTL.vtlInitialize(speaker_file_name)
failure = VTL.vtlSegmentSequenceToGesturalScore(segment_file_name, gesture_file_name)
if failure != 0:
    raise ValueError('Error in vtlSegmentSequenceToGesturalScore! Errorcode: %i' % failure)
failure = VTL.vtlGesturalScoreToAudio(gesture_file_name,  # input
                          wav_file_name,  # output can be empty string "", if
                                          # you don't want to save audio to file
                          None,  # *audio
                          None,  # *numSamples
                          1)  # enableConsoleOutput
if failure != 0:
    raise ValueError('Error in vtlGesturalScoreToAudio! Errorcode: %i' % failure)
VTL.vtlClose()


print("Part 2: export ema points only")
speaker_file_name = ctypes.c_char_p('../resources/JD3.speaker'.encode())
gesture_file_name = ctypes.c_char_p('apfelsine.ges'.encode())
ema_file_name = ctypes.c_char_p('apfelsine-ema.txt'.encode())
VTL.vtlInitialize(speaker_file_name)
failure = VTL.vtlGesturalScoreToEma(gesture_file_name,  # input
                          ema_file_name)  # output
if failure != 0:
    raise ValueError('Error in vtlGesturalScoreToEma! Errorcode: %i' % failure)
VTL.vtlClose()


print("Part 3: export mesh and ema points")
speaker_file_name = ctypes.c_char_p('../resources/JD3.speaker'.encode())
gesture_file_name = ctypes.c_char_p('apfelsine.ges'.encode())
file_name = ctypes.c_char_p('apfelsine'.encode())
path_name = ctypes.c_char_p('Meshes/apfelsine'.encode())
if os.path.exists(path_name.value.decode()):
    shutil.rmtree(path_name.value.decode())
os.mkdir(path_name.value.decode())
VTL.vtlInitialize(speaker_file_name)
failure = VTL.vtlGesturalScoreToEmaAndMesh(gesture_file_name,  # input
                          path_name, file_name)  # output
VTL.vtlClose()

if failure != 0:
    raise ValueError('Error in vtlGesturalScoreToEmaAndMesh! Errorcode: %i' % failure)


print("Part 4: Visualize mesh and ema points interactively")
print("NOTE: for this part you need `blender` and the python packages "
      "`pyglet` and `pywavefront`.")

# You need to have generated the meshes and ema export with Part 3!

# change the `convert_directory.py` to your liking and then execute it on the
# Terminal with blender
if sys.platform.startswith('win32'):
    os.system(r'"D:\Program Files\Blender Foundation\Blender 3.0\blender.exe" --background --python Meshes/convert_directory.py')
else:
    os.system('blender --background --python Meshes/convert_directory.py 1> /dev/null')

# now you can visualize the meshes interactively by executing the animation.py
if sys.platform.startswith('win32'):
    os.system('python.exe Meshes/animation.py')
else:
    os.system('python Meshes/animation.py')

# this pops up a new window with the animation and gives some output in the Terminal.

