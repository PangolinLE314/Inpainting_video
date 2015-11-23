% This script initialise inpainting at coarsest level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      INITIALISATION       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

occVolIter = occVol;
sumHole = sum(occVolIter(:)); %sum all the value in occlusion. (number of occlusion pixel)

X= sprintf('Initialisation for coarsest level\n--------------------------');
disp(X)
numOnionPeel=0;
%fill in, in an onion peel fashion
while(sumHole >0)              % fill until sum=0
    numOnionPeel=numOnionPeel+1;
    X= sprintf('Processing inpainting onion level %d \n--------------------------',numOnionPeel);
disp(X)
    sumHole = sum(occVolIter(:));
    
    structElCube = strel('arbitrary', ones(3,3,3)); %Create morphological structuring element, 8-square neighborhood
    occVolErode = imerode(occVolIter,structElCube); %erode image
    
    %set up the partial occlusion volume for the PatchMatch :
    % - 0 for non-occlusion;
    % - 1 for occluded and not to take into account when
    % comparing patches
    % - 2 for occluded and to take into account (we do not allow
    %the nearest neighbours to point to these pixels, but
    %they have been reconstructed at this iteration
    occVolPatchMatch = occVolDilate;
    occVolPatchMatch((occVolDilate - occVolIter) == 1) = 2;
    
    %initial guess : by default, set everything to 0
    if (~exist('shiftVol'))
        firstGuess = single(zeros([4 size(imgVol,2) size(imgVol,3) size(imgVol,4)]));
    else
        firstGuess = shiftVol+1-1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %carry out the 3D patchMatch
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate PatchMatch3D from imgVol
    shiftVol = spatio_temporal_patch_match_mex(imgVol, imgVol,...
        patchMatchParams,firstGuess,occVolPatchMatch,occVolDilate);
    
    if (exist('stop_and_debug.txt'))
        keyboard;
    end
    
    %occVolErode_show=occVolErode(:,:,1);
    %occVolBorder_show=occVolBorder(:,:,1)
    %determine the pixels to reconstruct at this layer
    occVolBorder = abs(occVolIter - occVolErode);
    %if the occVol == 2, then we cannot use the colour
    %info, but we do not inpaint it either
    occVolBorder(occVolErode == 1) = 2;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %carry out the reconstruction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (exist('featurePyramid','var') && exist('shiftVol'))
        [imgVol,gradX,gradY,normGradX,normGradY] = reconstruct_video_and_features_mex(imgVol,occVolBorder,...
            shiftVol,patchMatchParams,...
            sigmaColour,useAllPatches,reconstructionType);%
        patchMatchParams.gradX = single(gradX); patchMatchParams.gradY = single(gradY); patchMatchParams.normGradX = single(normGradX); patchMatchParams.normGradY = single(normGradY);
    else
        [imgVol] = reconstruct_video_mex(imgVol,occVolBorder,...
            shiftVol,patchMatchParams,...
            sigmaColour,useAllPatches,reconstructionType);%
    end
    
    %now go to the next onion peel layer
    occVolIter = occVolErode;
    %add 1 to iteration number
    iterationNb = iterationNb+1;
end

disp('Finish intialisation')
