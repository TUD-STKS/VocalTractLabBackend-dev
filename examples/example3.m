%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example simply shows how to obtain the volume velocity transfer
% function of the vocal tract based on vocal tract parameters for a certain
% phone in the speaker file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('***************** Running example3.m *****************')

if isfolder('../include/VocalTractLabApi')
    addpath('../include/VocalTractLabApi');
end
if isfolder('../lib/Release')
    addpath('../lib/Release');
end
%% Copy speaker file to examples folder
% This is necessary because folders named 'resources' cannot be added to
% the MATLAB path since R2019
status = copyfile('../resources/JD3.speaker', 'JD3.speaker');

%%
libName = 'VocalTractLabApi';
if isunix
    libName = ['lib' libName];
end
headerName = 'VocalTractLabApi.h';
    
if ~libisloaded(libName)
    % To load the library, specify the name of the DLL and the name of the
    % header file. If no file extensions are provided (as below)
    % LOADLIBRARY assumes that the DLL ends with .dll or .so
    loadlibrary(libName, headerName);
    disp(['Loaded library: ' libName]);
    pause(1);
end

if ~libisloaded(libName)
    error('Failed to load external library: %s', libName);
end

% *****************************************************************************
% list the methods
% *****************************************************************************

libfunctions(libName);   

% *****************************************************************************
% Print the version (compile date) of the library.
%
% void vtlGetVersion(char *version);
% *****************************************************************************

% Init the variable version with enough characters for the version string
% to fit in.
version = '                                ';
version = calllib(libName, 'vtlGetVersion', version);

disp(['Compile date of the library: ' version]);

% *****************************************************************************
% Initialize the VTL synthesis with the given speaker file name.
%
% void vtlInitialize(const char *speakerFileName)
% *****************************************************************************

speakerFileName = 'JD3.speaker';

failure = calllib(libName, 'vtlInitialize', speakerFileName);
if (failure ~= 0)
    error('Error in vtlInitialize()! Error code: %d', failure);   
end

% *****************************************************************************
% Get some constants.
% *****************************************************************************

audioSamplingRate = 0;
numTubeSections = 0;
numVocalTractParams = 0;
numGlottisParams = 0;
numAudioSamplesPerTractState = 0;
internalSamplingRate = 0;

[failure, audioSamplingRate, numTubeSections, numVocalTractParams, ...
    numGlottisParams, numAudioSamplesPerTractState, ...
    internalSamplingRate] = ...
    calllib(libName, 'vtlGetConstants', ...
    audioSamplingRate, numTubeSections, numVocalTractParams, ...
    numGlottisParams, numAudioSamplesPerTractState, internalSamplingRate);

if (failure ~= 0)
    error('Error in vtlGetConstants()! Error code: %d', failure);   
end

% *****************************************************************************
% Get the vocal tract parameters for the phone /i/.
%
% int vtlGetTractParams(char *shapeName, double *param);
% *****************************************************************************

vocalTractParams = zeros(1, numVocalTractParams);
shapeName = 'i';

[failure, shapeName, vocalTractParams] = ...
  calllib(libName, 'vtlGetTractParams', shapeName, vocalTractParams);

if (failure ~= 0)
    error('Error in vtlGetTractParams()! Error code: %d', failure);   
end


%%
opts.spectrumType = 'SPECTRUM_UU';
opts.radiationType = 'NO_RADIATION';
opts.boundaryLayer = false;
opts.heatConduction = false;
opts.softWalls = false;
opts.hagenResistance = false;
opts.innerLengthCorrections = false;
opts.lumpedElements = false;
opts.paranasalSinuses = false;
opts.piriformFossa = false;
opts.staticPressureDrops = false;

[failure, opts] = ...
    calllib(libName, 'vtlGetDefaultTransferFunctionOptions', ...
    opts);

if (failure ~= 0)
    error('Error in vtlGetDefaultTransferFunctionOptions()! Error code: %d', failure);   
end


%%
% *****************************************************************************
% int vtlGetTransferFunction(double* tractParams, int numSpectrumSamples,
%    TransferFunctionOptions* opts, double* magnitude, double* phase_rad);
% *****************************************************************************
NUM_SPECTRUM_SAMPLES = 2048;
magSpectrum = zeros(1, NUM_SPECTRUM_SAMPLES);
phaseSpectrum = zeros(1, NUM_SPECTRUM_SAMPLES);

[failure, vocalTractParams, opts, magSpectrum, phaseSpectrum] = ...
  calllib(libName, 'vtlGetTransferFunction', vocalTractParams, ...
    NUM_SPECTRUM_SAMPLES, opts, magSpectrum, phaseSpectrum);

if (failure ~= 0)
    error('Error in vtlGetTransferFunction()! Error code: %d', failure);   
end

% Plot the transfer function up to 10000 Hz.

numPlotSamples = int32(10000 * NUM_SPECTRUM_SAMPLES / audioSamplingRate);
freqAxis = double(0:1:numPlotSamples-1);
freqAxis = (double(audioSamplingRate) / double(NUM_SPECTRUM_SAMPLES)).*freqAxis;
plot(freqAxis, 20*log10(magSpectrum(1:numPlotSamples)), freqAxis, phaseSpectrum(1:numPlotSamples));

%ylim([-40, 40]);
xlabel('Frequency in Hz');
ylabel('Log. magnitude in dB, phase in rad');

%%
% *****************************************************************************
% Close the VTL synthesis.
%
% void vtlClose();
% *****************************************************************************

calllib(libName, 'vtlClose');
if (failure ~= 0)
    error('Error in vtlClose()! Error code: %d', failure);
end

unloadlibrary(libName);

disp('***************** Finished example3.m *****************')
