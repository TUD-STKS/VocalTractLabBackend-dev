#!/usr/bin/env python3

'''
This example generates the speech waveform directly from a gestural score.

Look into example1.py for more thorough comments on how to interface
vocaltractlab API from python3.

'''

import ctypes
import sys

# Use 'VocalTractLabApi32.dll' if you use a 32-bit python version.

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

# Synthesize from a gestural score.
gesture_file_name = ctypes.c_char_p(b'ala.ges')
wav_file_name = ctypes.c_char_p(b'ala.wav')
num_samples = ctypes.c_int(0)
audio = (ctypes.c_double * int(1*44100))()  # reserve memory

print('Calling vtlGesturalScoreToAudio()...')

failure = VTL.vtlGesturalScoreToAudio(gesture_file_name,
                                      wav_file_name,  # can be "" to not save anything
                                      ctypes.byref(audio),  # can be None to not store audio
                                      ctypes.byref(num_samples),
                                      1)  # can be 0 to suppress output

if failure != 0:
    raise ValueError('Error in vtlGesturalScoreToAudio! Errorcode: %i' % failure)

failure = VTL.vtlClose()

print('Finished.')

