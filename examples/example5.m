%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example simply shows how to 
% o transform a segment sequence file into a gestural score file,
% o transform the gestural score file into a tract parameter sequence file,
% o generate audio from the tract parameter sequence file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('***************** Running example5.m *****************')

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

%%
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
% int vtlSegmentSequenceToGesturalScore(const char *segFileName, ...
%   const char *gesFileName, bool enableConsoleOutput);
% *****************************************************************************

segFileName = 'apfelsine.seg';
gesFileName = 'apfelsine.ges';
enableConsoleOutput = true;
tractSequenceFileName = 'apfelsine.txt';
wavFileName = 'apfelsine.wav';

failure = ...
  calllib(libName, 'vtlSegmentSequenceToGesturalScore', ...
  segFileName, gesFileName, enableConsoleOutput);

if (failure ~= 0)
    error('Error in vtlSegmentSequenceToGesturalScore()! Error code: %d', failure);   
end

audio = zeros(10*44100, 1);     % Reserve memory for 10 s
numSamples = 0;

% *****************************************************************************
% int vtlGesturalScoreToAudio(const char* gesFileName, const char* wavFileName,
%  double* audio, int* numSamples, int enableConsoleOutput);
% *****************************************************************************

% Possibly synthesize the audio directly from the gestural score...

%[failure, gesFileName, wavFileName, audio, numSamples] = ...
%  calllib(libName, 'vtlGesturalScoreToAudio', gesFileName, wavFileName, audio, numSamples, 1);
%if (failure ~= 0)
%    error('Error in vtlGesturalScoreToAudio()! Error code: %d', failure);   
%end
%plot(1:numSamples, audio(1:numSamples));

% *****************************************************************************
% int vtlGesturalScoreToTractSequence(const char* gesFileName, 
%  const char* tractSequenceFileName);
% *****************************************************************************

failure = ...
  calllib(libName, 'vtlGesturalScoreToTractSequence', gesFileName, tractSequenceFileName);

if (failure ~= 0)
    error('Error in vtlGesturalScoreToTractSequence()! Error code: %d', failure);   
end

% *****************************************************************************
% *****************************************************************************

% Here, the sequence of parameter values of the vocal tract or glottis models 
% could be manipulated in the tract sequence TXT file.

% *****************************************************************************
% int vtlTractSequenceToAudio(const char* tractSequenceFileName,
%  const char* wavFileName, double* audio, int* numSamples);
% *****************************************************************************

[failure, tractSequenceFileName, wavFileName, audio, numSamples] = ...
  calllib(libName, 'vtlTractSequenceToAudio', tractSequenceFileName, wavFileName, audio, numSamples);

if (failure ~= 0)
    error('Error in vtlTractSequenceToAudio()! Error code: %d', failure);   
end

plot(1:numSamples, audio(1:numSamples));


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

disp('***************** Finished example5.m *****************')