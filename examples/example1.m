%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example generates the transition from /a/ to /i/ using the 
% incremental tract synthesis with the function 
% vtlSynthesisAddTract(...).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('***************** Running example1.m *****************')

%% Add paths containing header file and binary of library
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

%%
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

%%
% *****************************************************************************
% Initialize the VTL synthesis with the given speaker file name.
%
% void vtlInitialize(const char *speakerFileName)
% *****************************************************************************

speakerFileName = 'JD3.speaker';

failure = calllib(libName, 'vtlInitialize', which(speakerFileName));
if (failure ~= 0)
    error('Error in vtlInitialize()! Error code: %d', failure);   
end

%%
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

disp(['Audio sampling rate = ' num2str(audioSamplingRate)]);
disp(['Num. of tube sections = ' num2str(numTubeSections)]);
disp(['Num. of vocal tract parameters = ' num2str(numVocalTractParams)]);
disp(['Num. of glottis parameters = ' num2str(numGlottisParams)]);
disp(['Num. of audio samples per tract state = ' ...
    num2str(numAudioSamplesPerTractState)]);
disp(['Internal sampling rate = ', num2str(internalSamplingRate)]);

%%
% *****************************************************************************
% Get information about the parameters of the vocal tract model and the
% glottis model.
%
% void vtlGetTractParamInfo(char *names, double *paramMin, double *paramMax, 
%   double *paramNeutral);
% void vtlGetGlottisParamInfo(char *names, double *paramMin, double *paramMax, 
%   double *paramNeutral);
% *****************************************************************************

% Reserve enough chars for each parameter.
tractParamNames = blanks(numVocalTractParams*10);
tractParamDescriptions = blanks(numVocalTractParams*100);
tractParamUnits = blanks(numVocalTractParams*10);
tractParamMin = zeros(1, numVocalTractParams);
tractParamMax = zeros(1, numVocalTractParams);
tractParamNeutral = zeros(1, numVocalTractParams);

[failure, tractParamNames, tractParamDescriptions, tractParamUnits, ...
    tractParamMin, tractParamMax, tractParamNeutral] = ...
  calllib(libName, 'vtlGetTractParamInfo', ...
  tractParamNames, tractParamDescriptions, tractParamUnits, ...
  tractParamMin, tractParamMax, tractParamNeutral);

if (failure ~= 0)
    error('Error in vtlGetTractParamInfo()! Error code: %d', failure);   
end

disp(['Vocal tract parameter names: ' tractParamNames]);
disp(['Vocal tract parameter descriptions: ' tractParamDescriptions]);

%%
glottisParamNames = blanks(numGlottisParams*10);
glottisParamDescriptions = blanks(numGlottisParams*100);
glottisParamUnits = blanks(numGlottisParams*10);
glottisParamMin = zeros(1, numGlottisParams);
glottisParamMax = zeros(1, numGlottisParams);
glottisParamNeutral = zeros(1, numGlottisParams);

[failure, glottisParamNames, glottisParamDescriptions, glottisParamUnits, ...
    glottisParamMin, glottisParamMax, glottisParamNeutral] = ...
  calllib(libName, 'vtlGetGlottisParamInfo', ...
  glottisParamNames, glottisParamDescriptions, glottisParamUnits, ...
  glottisParamMin, glottisParamMax, glottisParamNeutral);

if (failure ~= 0)
    error('Error in vtlGetGlottisParamInfo()! Error code: %d', failure);   
end

disp(['Glottis parameter names: ' glottisParamNames]);
disp(['Glottis parameter descriptions: ' glottisParamDescriptions]);
%%
% *****************************************************************************
% Get the vocal tract parameter values for the vocal tract shapes of /i/
% and /a/, which are saved in the speaker file.
%
% int vtlGetTractParams(char *shapeName, double *param);
% *****************************************************************************

shapeName = 'a';
paramsA = zeros(1, numVocalTractParams);
[failed, ~, paramsA] = ...
  calllib(libName, 'vtlGetTractParams', shapeName, paramsA);

if (failed ~= 0)
    error('Error: Vocal tract shape "%s" not in the speaker file!', shapeName);   
end

shapeName = 'i';
paramsI = zeros(1, numVocalTractParams);
[failed, shapeName, paramsI] = ...
  calllib(libName, 'vtlGetTractParams', shapeName, paramsI);

if (failed ~= 0)
    error('Error: Vocal tract shape "%s" not in the speaker file!', shapeName);   
end
%%
% *****************************************************************************
% Incrementally synthesize a transition from /a/ to /i/.
%
% void vtlSynthesisReset()
%
% int vtlSynthesisAddTract(int numNewSamples, double *audio,
%   double *tractParams, double *glottisParams);
% *****************************************************************************

audio1 = zeros(1, 10000);
audio2 = zeros(1, 10000);
audio3 = zeros(1, 10000);
audio4 = zeros(1, 10000);

glottisParams = glottisParamNeutral;

% Initialize the tube synthesis.
calllib(libName, 'vtlSynthesisReset');

% Submit the initial vocal tract shape (numSamples=0) with P_sub = 0
glottisParams(2) = 0;       % P_sub = 0 dPa
[failure, audio1, ~, glottisParams] = ...
  calllib(libName, 'vtlSynthesisAddTract', 0, audio1, ...
    paramsI, glottisParams);

if (failure ~= 0)
    error('Error in vtlSynthesisAddTract()! Error code: %d', failure);   
end

% Ramp up the subglottal pressure within 1000 samples
glottisParams(2) = 8000;   % P_sub = 8000 dPa
[failure, audio1, ~, glottisParams] = ...
  calllib(libName, 'vtlSynthesisAddTract', 1000, audio1, ...
    paramsI, glottisParams);

if (failure ~= 0)
    error('Error in vtlSynthesisAddTract()! Error code: %d', failure);   
end

% Make transitions between /a/ and /i/
[failure, audio2, ~, glottisParams] = ...
  calllib(libName, 'vtlSynthesisAddTract', 10000, audio2, ...
    paramsA, glottisParams);

if (failure ~= 0)
    error('Error in vtlSynthesisAddTract()! Error code: %d', failure);   
end

[failure, audio3, ~, glottisParams] = ...
  calllib(libName, 'vtlSynthesisAddTract', 10000, audio3, ...
    paramsI, glottisParams);

if (failure ~= 0)
    error('Error in vtlSynthesisAddTract()! Error code: %d', failure);   
end

[failure, audio4, tractParams, glottisParams] = ...
  calllib(libName, 'vtlSynthesisAddTract', 10000, audio4, ...
    paramsA, glottisParams);

if (failure ~= 0)
    error('Error in vtlSynthesisAddTract()! Error code: %d', failure);   
end

audio = [audio1(1:1000) audio2 audio3 audio4];

% Plot and play the audio signal
plot(audio);
soundsc(audio, double(audioSamplingRate));
audiowrite('test.wav', audio, audioSamplingRate);


% *****************************************************************************
% Close the VTL synthesis.
%
% void vtlClose();
% *****************************************************************************
disp('Calling vtlClose...');
calllib(libName, 'vtlClose');
if (failure ~= 0)
    error('Error in vtlClose()! Error code: %d', failure);
end

unloadlibrary(libName);
disp('***************** Finished example1.m *****************')
