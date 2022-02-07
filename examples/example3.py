#!/usr/bin/env python3

'''
This example shows how to obtain the volume velocity transfer function of the
vocal tract based on vocal tract parameters for a certain phone in the speaker
file.

Look into example1.py for more thorough comments on how to interface
vocaltractlab API from python3.

'''

import ctypes
import sys

# try to load some non-essential packages
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
try:
    import numpy as np
except ImportError:
    np = None


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


# get some constants
audio_sampling_rate = ctypes.c_int(0)
number_tube_sections = ctypes.c_int(0)
number_vocal_tract_parameters = ctypes.c_int(0)
number_glottis_parameters = ctypes.c_int(0)
number_audio_samples_per_tract_state = ctypes.c_int(0)
internal_sampling_rate = ctypes.c_double(0)

failure = VTL.vtlGetConstants(
        ctypes.byref(audio_sampling_rate),
        ctypes.byref(number_tube_sections),
        ctypes.byref(number_vocal_tract_parameters),
        ctypes.byref(number_glottis_parameters),
        ctypes.byref(number_audio_samples_per_tract_state),
        ctypes.byref(internal_sampling_rate))

if failure != 0:
    raise ValueError('Error in vtlGetConstants! Errorcode: %i' % failure)

print('Audio sampling rate = %i' % audio_sampling_rate.value)
print('Num. of tube sections = %i' % number_tube_sections.value)
print('Num. of vocal tract parameters = %i' % number_vocal_tract_parameters.value)
print('Num. of glottis parameters = %i' % number_glottis_parameters.value)
print('Num. of audio samples per tract state = %i' % number_audio_samples_per_tract_state.value)
print('Internal sampling rate = %d' % internal_sampling_rate.value)



# Get the vocal tract parameters for the phone /a/, which are saved in the
# speaker file.
TRACT_PARAM_TYPE = ctypes.c_double * number_vocal_tract_parameters.value
shape_name = ctypes.c_char_p(b'a')
params_a = TRACT_PARAM_TYPE()

failure = VTL.vtlGetTractParams(shape_name, ctypes.byref(params_a))
if failure != 0:
    raise ValueError('Error in vtlGetTractParams! Errorcode: %i' % failure)


# extract transfer function
NUM_SPECTRUM_SAMPLES = 1024
SPECTRUM_TYPE = ctypes.c_double * NUM_SPECTRUM_SAMPLES
magnitude_spectrum = SPECTRUM_TYPE()
phase_spectrum = SPECTRUM_TYPE()  # in radiants

VTL.vtlGetTransferFunction(ctypes.byref(params_a),  # input
                           NUM_SPECTRUM_SAMPLES,  # input
                           None,  # input TransferFunctionOptions
                           ctypes.byref(magnitude_spectrum),  # output
                           ctypes.byref(phase_spectrum))  # output

print('First 40 data points for every vector:')
print('  magnitude_spectrum: %s' % str(list(magnitude_spectrum)[:40]))
print('  phase_spectrum: %s' % str(list(phase_spectrum)[:40]))


# destroy current state of VTL and free memory
VTL.vtlClose()


# plot and play transfer function
#################################

if plt is not None and np is not None:
    frequency = np.arange(NUM_SPECTRUM_SAMPLES)
    frequency = frequency * (float(audio_sampling_rate.value) / NUM_SPECTRUM_SAMPLES)
    plt.plot(frequency, 20 * np.log10(magnitude_spectrum), c='black',
             linewidth=2.0, alpha=0.75, label='magnitude [dB]')
    plt.plot(frequency, phase_spectrum, c='red', linewidth=2.0,
             linestyle=':', alpha=0.75, label='phase [rad]')
    plt.ylabel('magnitude [dB], phase [rad]')
    plt.xlabel('frequency [Hz]')
    plt.xlim(0,10000)
    plt.ylim(-50,50)
    plt.title('volume velocity transfer function')
    plt.legend()
    plt.tight_layout()
    plt.show(block=False)
else:
    print('plotting not available; matplotlib needed')
    print('skip plotting')

