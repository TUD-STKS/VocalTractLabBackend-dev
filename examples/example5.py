#!/usr/bin/env python3

'''
This example simply shows how to

 - transform a segment sequence file into a gestural score file,
 - transform the gestural score file into a tract parameter sequence file,
 - generate audio from the tract parameter sequence file.

Look into example1.py for more thorough comments on how to interface
vocaltractlab API from python3.

'''

import ctypes
import sys


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


# Synthesize from segment sequence via gestural score.
speaker_file_name = ctypes.c_char_p(b'../resources/JD3.speaker')
segment_file_name = ctypes.c_char_p(b'apfelsine.seg')
gesture_file_name = ctypes.c_char_p(b'apfelsine.ges')
tract_sequence_file_name = ctypes.c_char_p(b'apfelsine.txt')
wav_file_name = ctypes.c_char_p(b'apfelsine.wav')

print("VTL.vtlInitialize...")
failure = VTL.vtlInitialize(speaker_file_name)
if failure != 0:
    raise ValueError('Error in vtlInitialize! Errorcode: %i' % failure)

print("VTL.vtlSegmentSequenceToGesturalScore...")
failure = VTL.vtlSegmentSequenceToGesturalScore(segment_file_name, gesture_file_name)
if failure != 0:
    raise ValueError('Error in vtlSegmentSequenceToGesturalScore! Errorcode: %i' % failure)

print("VTL.vtlGesturalScoreToTractSequence...")
failure = VTL.vtlGesturalScoreToTractSequence(gesture_file_name, tract_sequence_file_name)
if failure != 0:
    raise ValueError('Error in vtlGesturalScoreToTractSequence! Errorcode: %i' % failure)

print("VTL.vtlTractSequenceToAudio...")
failure = VTL.vtlTractSequenceToAudio(tract_sequence_file_name, wav_file_name, None, None)
if failure != 0:
    raise ValueError('Error in vtlTractSequenceToAudio! Errorcode: %i' % failure)

failure = VTL.vtlClose()
if failure != 0:
    raise ValueError('Error in vtlClose! Errorcode: %i' % failure)

print('Finished.')

