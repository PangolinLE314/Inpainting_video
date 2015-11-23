%this script creat inpainting video script for debug.
% Date : 06/11/2015
% Author: LE Thuc Trinh.
close all;
clear all;
name='truck';
restoredefaultpath;
% videoFile = './Content/ball.avi';
% occlusionFile = './Content/rubikOcc.avi';
% outputVideoFile='rubik_inpainted.avi';

videoFile = strcat('./Content/',name,'.avi');
occlusionFile = strcat('./Content/',name,'Occ.avi');
outputVideoFile=strcat(name,'_inpainted.avi');

[videoPath,file,~] = fileparts(videoFile);
[occlusionPath,~,~] = fileparts(occlusionFile);
currDir = cd;
addpath(genpath(currDir));
%turn off warning if it exists
warning('off','MATLAB:aviread:FunctionToBeRemoved');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% get the input video %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input video : imgVol 4-D matrix: 1,2: Image, 3:3 color channel , 4:Time.
disp('Reading input video');
imgVol = read_video(videoFile);
if (isempty(imgVol))
    return;
end
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










