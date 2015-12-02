%this script creat inpainting video script for debug.
% Date : 06/11/2015
% Author: LE Thuc Trinh.
%close all;
%clear all;
name='chessboard';
rootDataset='D:\PhD\Results\Videos\';
restoredefaultpath;
videoFile = fullfile(rootDataset,strcat(name,'.avi'));
occlusionFile = fullfile(rootDataset,strcat(name,'Occ.avi'));
outputVideoFile=strcat(name,'_inpainted.avi');
%turn off warning if it exists
warning('off','MATLAB:aviread:FunctionToBeRemoved');

[videoPath,file,~] = fileparts(videoFile);
[occlusionPath,~,~] = fileparts(occlusionFile);
currDir = cd;
addpath(genpath(currDir));
addpath(videoPath,occlusionPath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% get the input video %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input video : imgVol 4-D matrix: 1,2: Image, 3:3 color channel , 4:Time.
disp('Reading input video');
imgVol = read_video(videoFile);
if (isempty(imgVol))
    return;
end
imgVol=uint8(imgVol);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% get the occulusion video %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Reading input occlusion');
occVol = read_video(occlusionFile);
%case of colour occlusion
if (ndims(occVol) == 4)
    occVol(:,:,:,2:3) = [];
end
% special case (quite common) where the occlusion file is one image
if (ndims(occVol) == 2)
    occVol = repmat(occVol,[1 1 size(imgVol,3)]);
    if (isempty(occVol))
        return;
    end
end
occVol= uint8(occVol);
%check the dimensions of the videos
sizeVideo = size(imgVol);
sizeOccVol =size(occVol);
if ( ~(isequal(sizeVideo(1:2),sizeOccVol(1:2))))
    disp('problem, the input video and occlusion do not have the same dimension. exiting the program.');
    return;
end
if (size(imgVol,3) > 3)	%we need to permute the dimensions
    imgVol = permute(imgVol,[1 2 4 3]);
end

if ( size(imgVol,4) ~= size(occVol,3) )
    disp('There was an error reading the input : the number of frames in the video and the occlusion are not the same');
    return;
end

%% Initialise params for inpainting
disp('Starting the video inpainting ! Hold on to your hats !!!');
profile clear;
profile on;

parameter = {videoFile,occlusionFile,'maxLevel',3,'textureFeaturesActivated',1,'patchSizeX',ps,'patchSizeY',ps,'patchSizeT',ps};
inpainting_parameters = parameter(3:end);
inpainting_parameters{end+1} = 'file';
inpainting_parameters{end+1} = file;
imgVolIn= single(permute(imgVol,[3 1 2 4]));
occVolIn= occVol;

%FIXED PARMETERS !!!
%GAUSSIAN PYRAMID PARAMETERS
filterSize = 3;    %fixed
scaleStep = 0.5;    %fixed
sigma = 1.5;    %fixed
useAllPatches = 0;
reconstructionType = 0;
[maxLevel,patchSizeX,patchSizeY,patchSizeT,textureFeaturesActivated,sigmaColour,file] = ...
    parse_inpaint_parameters(inpainting_parameters);
%patchMatch params
patchSize = [patchSizeX patchSizeY min(patchSizeT,size(imgVolIn,4))]; %size of small cube to search in PatchMatch
patchMatchParams.patchSize = [patchSizeX patchSizeY min(patchSizeT,size(imgVolIn,4))];
patchMatchParams.patchSizeX = patchSize(1);
patchMatchParams.patchSizeY = patchSize(2);
patchMatchParams.patchSizeT = patchSize(3);
patchMatchParams.w = max([size(imgVolIn,2) size(imgVolIn,3) size(imgVolIn,4)]);  % maximum width (not sure)
patchMatchParams.alpha = 0.5;   %fixed
patchMatchParams.fullSearch = 0;
patchMatchParams.partialComparison = 1;
patchMatchParams.nbItersPatchMatch = 10;
patchMatchParams.patchIndexing = 0;
patchMatchParams.reconstructionType = reconstructionType;

%parameters concerning the iterations of the nearest neighbour search
%and reconstruction
maxNbIterations = 2;
residualThresh = 0.1;

%%
script_inpaiting_video;









