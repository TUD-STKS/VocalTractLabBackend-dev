%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example generates the speech waveform directly from a gestural
% score.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('***************** Running example2.m *****************')

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
% Synthesize from a gestural score.
%
% int vtlGesturalScoreToWav(const char *gesFileName, const char *wavFileName,
%  int enableConsoleOutput);
% *****************************************************************************

speakerFileName = 'JD3.speaker';
gestureFileName = 'ala.ges';
wavFileName = 'ala.wav';

failure = calllib(libName, 'vtlInitialize', speakerFileName);

if (failure ~= 0)
    error('Error in vtlInitialize()! Error code: %d', failure);
end

numSamples = 0;
audio = zeros(44100, 0);   % Enough for 1 s of audio.

failure = calllib(libName, 'vtlGesturalScoreToAudio', gestureFileName, ...
    wavFileName, audio, numSamples, 1);

if (failure ~= 0)
    error('Error in vtlGesturalScoreToAudio()! Error code: %d', failure);
end

failure = calllib(libName, 'vtlClose');
if (failure ~= 0)
    error('Error in vtlClose()! Error code: %d', failure);
end
unloadlibrary(libName);
disp('Finished.');

% Play the synthesized wav file.
s = audioread(wavFileName);

plot(1:length(s), s);
sound(s, 44100);

disp('***************** Finished example2.m *****************')

