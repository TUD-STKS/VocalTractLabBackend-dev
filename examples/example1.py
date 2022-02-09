#!/usr/bin/env python3

'''
This example generates the transition from /a/ to /i/ using the vocal tract
model and the function vtlSynthesisAddTract(...).

.. note::

    This example uses ``ctypes`` and is very close to the vocaltractlab API.
    This comes with breaking with some standard assumptions one have in python.

If you are not aware of ``ctypes`` read the following introduction
https://docs.python.org/3/library/ctypes.html

For an in-depth API description look at the `../include/VocalTractLabApi/VocalTractLabApi.h`.

For plotting and saving results you need to install ``matplotlib``, ``numpy``,
and ``scipy``.

'''

import ctypes
import sys

# try to load some non-essential packages
try:
    import numpy as np
except ImportError:
    np = None
try:
    from scipy.io import wavfile
except ImportError:
    wavfile = None
try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


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


# get information about the parameters of the vocal tract model
# Hint: Reserve 10; 100; 10 chars for each parameter.
TRACT_PARAM_TYPE = ctypes.c_double * number_vocal_tract_parameters.value
tract_param_names = ctypes.c_char_p((' ' * 10 * number_vocal_tract_parameters.value).encode())
tract_param_descriptions = ctypes.c_char_p((' ' * 100 * number_vocal_tract_parameters.value).encode())
tract_param_units = ctypes.c_char_p((' ' * 10 * number_vocal_tract_parameters.value).encode())
tract_param_min = TRACT_PARAM_TYPE()
tract_param_max = TRACT_PARAM_TYPE()
tract_param_neutral = TRACT_PARAM_TYPE()

VTL.vtlGetTractParamInfo(tract_param_names,
                         tract_param_descriptions,
                         tract_param_units,
                         ctypes.byref(tract_param_min),
                         ctypes.byref(tract_param_max),
                         ctypes.byref(tract_param_neutral))

print('Vocal tract parameters: "%s"' % tract_param_names.value.decode())
print('Vocal tract descriptions: "%s"' % tract_param_descriptions.value.decode())
print('Vocal tract units: "%s"' % tract_param_units.value.decode())
print('Vocal tract parameter minima: ' + str(list(tract_param_min)))
print('Vocal tract parameter maxima: ' + str(list(tract_param_max)))
print('Vocal tract parameter neutral: ' + str(list(tract_param_neutral)))

# get information about the parameters of glottis model
# Hint: Reserve 10; 100; 10 chars for each parameter.
GLOTTIS_PARAM_TYPE = ctypes.c_double * number_glottis_parameters.value
glottis_param_names = ctypes.c_char_p((' ' * 10 * number_glottis_parameters.value).encode())
glottis_param_descriptions = ctypes.c_char_p((' ' * 100 * number_glottis_parameters.value).encode())
glottis_param_units = ctypes.c_char_p((' ' * 10 * number_glottis_parameters.value).encode())
glottis_param_min = GLOTTIS_PARAM_TYPE()
glottis_param_max = GLOTTIS_PARAM_TYPE()
glottis_param_neutral = GLOTTIS_PARAM_TYPE()

VTL.vtlGetGlottisParamInfo(glottis_param_names,
                           glottis_param_descriptions,
                           glottis_param_units,
                           ctypes.byref(glottis_param_min),
                           ctypes.byref(glottis_param_max),
                           ctypes.byref(glottis_param_neutral))

print('Glottis parameters: "%s"' % glottis_param_names.value.decode())
print('Glottis descriptions: "%s"' % glottis_param_descriptions.value.decode())
print('Glottis units: "%s"' % glottis_param_units.value.decode())
print('Glottis parameter minima: ' + str(list(glottis_param_min)))
print('Glottis parameter maxima: ' + str(list(glottis_param_max)))
print('Glottis parameter neutral: ' + str(list(glottis_param_neutral)))


# Get the vocal tract parameter values for the vocal tract shapes of /i/
# and /a/, which are saved in the speaker file.

shape_name = ctypes.c_char_p(b'a')
params_a = TRACT_PARAM_TYPE()
failure = VTL.vtlGetTractParams(shape_name, ctypes.byref(params_a))
if failure != 0:
    raise ValueError('Error in vtlGetTractParams! Errorcode: %i' % failure)

shape_name = ctypes.c_char_p(b'i')
params_i = TRACT_PARAM_TYPE()
failure = VTL.vtlGetTractParams(shape_name, ctypes.byref(params_i))
if failure != 0:
    raise ValueError('Error in vtlGetTractParams! Errorcode: %i' % failure)

del shape_name


# *****************************************************************************
# Incrementally synthesize a transition from /i/ to /a/ to /i/.
#
# void vtlSynthesisReset()
#
# int vtlSynthesisAddTract(int numNewSamples, double *audio,
#   double *tractParams, double *glottis_params);
#
#  example4.py shows how to do it with vtlSynthBlock in one sweep.
#
# *****************************************************************************

audios = [(ctypes.c_double * int(1000))(),
          (ctypes.c_double * int(10000))(),
          (ctypes.c_double * int(10000))(),
          (ctypes.c_double * int(10000))()]


glottis_params = glottis_param_neutral


# Initialize the tube synthesis.
VTL.vtlSynthesisReset()


# Submit the initial vocal tract shape (numSamples=0) with P_sub = 0
glottis_params[1] = 0  # P_sub = 0 dPa

VTL.vtlSynthesisAddTract(0, ctypes.byref(audios[0]), ctypes.byref(params_i), ctypes.byref(glottis_params))


# Ramp up the subglottal pressure within 1000 samples
glottis_params[1] = 8000  # P_sub = 8000 dPa
VTL.vtlSynthesisAddTract(1000, ctypes.byref(audios[0]), ctypes.byref(params_i), ctypes.byref(glottis_params))

# Make transitions between /a/ and /i/
VTL.vtlSynthesisAddTract(10000, ctypes.byref(audios[1]), ctypes.byref(params_a), ctypes.byref(glottis_params))

VTL.vtlSynthesisAddTract(10000, ctypes.byref(audios[2]), ctypes.byref(params_i), ctypes.byref(glottis_params))

VTL.vtlSynthesisAddTract(10000, ctypes.byref(audios[3]), ctypes.byref(params_a), ctypes.byref(glottis_params))

_wav = []
for wave in audios:
    _wav.extend(list(wave))


# destroy current state of VTL and free memory
VTL.vtlClose()


# plot and save the audio signal
################################

if np is not None:
    wav = np.array(_wav)

if plt is not None and np is not None:
    #time_wav = np.arange(len(wav)) / frame_rate
    plt.plot( wav, c='black', alpha=0.75)
    plt.ylabel('amplitude')
    plt.ylim(-1, 1)
    plt.xlabel('sample [s]')
    plt.show(block=False)
else:
    print('plotting not available; matplotlib needed')
    print('skip plotting')


if wavfile is not None and np is not None:
    wav_int = np.int16(wav  * (2 ** 15 - 1))
    wavfile.write('ai_test.wav', audio_sampling_rate.value, wav_int)
    print('saved audio to "ai_test.wav"')
else:
    print('scipy not available')
    print('skip writing out wav file')

# Plot the area function of the first and the last frame.
# TODO
# Matlab code:
#figure;
#plot(1:1:numTubeSections, tubeAreas(1:numTubeSections), ...
#    1:1:numTubeSections,
#    tubeAreas(1+(numFrames-1)*numTubeSections:(numFrames-1)*numTubeSections +
#    numTubeSections));
#xlabel('Position in cm');
#ylabel('Tube section index');

