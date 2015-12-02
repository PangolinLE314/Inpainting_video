%%
t1 = tic;

%%% get pyramid volumes
%imgVolPyramid: 3x3 cells store 3 chanels rgb.
imgVolPyramid = get_image_volume_pyramid(imgVolIn,filterSize,sigma,maxLevel,size(imgVolIn,4));%patchSize(3)
occVolPyramid = get_image_volume_pyramid(occVolIn,filterSize,sigma,maxLevel,size(imgVolIn,4));%patchSize(3)

%%%%% create texture feature pyramid
%descriptor using is causels et al, gradX,Y et normGradX,Y
if (textureFeaturesActivated>0)
    disp('Calculating texture feature pyramids');
    featurePyramid = get_video_features(imgVolIn,occVolIn,maxLevel,file);
end

clear imgVol;
%Start the loop:
for ii=maxLevel:-1:1                       %start with coarsest level
    X=sprintf('Processing level %d',ii); disp(X)
    pp=1;
    iterationNb = 1;
    residual = inf;
    
    occVol = single(occVolPyramid{ii});     %get the level of the pyramid
    
    imgVol(1,:,:,:) = imgVolPyramid{ii,1};  %get the red channel for this image volume
    imgVol(2,:,:,:) = imgVolPyramid{ii,2};  %get the blue channel for this image volume
    imgVol(3,:,:,:) = imgVolPyramid{ii,3};  %get the green channel for this image volume
    
    if (exist('featurePyramid','var'))
        gradX = single(featurePyramid{ii,1});
        gradY = single(featurePyramid{ii,2});
        normGradX = single(featurePyramid{ii,3});
        normGradY = single(featurePyramid{ii,4});
        
        patchMatchParams.gradX = gradX;
        patchMatchParams.gradY = gradY;
        patchMatchParams.normGradX = normGradX;
        patchMatchParams.normGradY = normGradY;
    end
    %if we are not at the coarsest level, then we recreate the image volume
    %using the previously upsampled shift map
    if (ii~=maxLevel)
        occVolToInpaint = occVol;
        if (exist('featurePyramid','var') && exist('shiftVol'))
            [imgVol,gradX,gradY,normGradX,normGradY] = reconstruct_video_and_features_mex(imgVol,occVolToInpaint,...
                shiftVol,patchMatchParams,...
                sigmaColour,useAllPatches,reconstructionType);%
            patchMatchParams.gradX = single(gradX); patchMatchParams.gradY = single(gradY); patchMatchParams.normGradX = single(normGradX); patchMatchParams.normGradY = single(normGradY);
        else
            [imgVol] = reconstruct_video_mex(imgVol,occVolToInpaint,...
                shiftVol,patchMatchParams,...
                sigmaColour,useAllPatches,reconstructionType);%
        end
    end
    
    %get the finer version of the image volume
    imgVolFine(1,:,:,:) = imgVolPyramid{max(ii-1,1),1};  %get the red channel for this image volume
    imgVolFine(2,:,:,:) = imgVolPyramid{max(ii-1,1),2};  %get the blue channel for this image volume
    imgVolFine(3,:,:,:) = imgVolPyramid{max(ii-1,1),3};  %get the green channel for this image volume
    
    %create a structuring element which is the size of the patch : this
    %will be used to dilate the occlusion
    structElPatch = strel('arbitrary', ones(patchSize(2),patchSize(1),patchSize(3)));
    occVolDilate = imdilate(occVol,structElPatch);
    
    iterationNb = iterationNb+1;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start of iterative inpainting at this level
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while (pp<= maxNbIterations && residual > residualThresh) % loop with maxIterations et residual < resThresh
        %%
        X=sprintf('=================\n Looping until convergence : %d',pp); disp(X)
        sizeImgVol = size(imgVol);
        imgVolIterMinusOne = imgVol+1-1; % question: what is this for?
        
        if (ii ~= maxLevel || pp>1)     %not bottom level
            patchMatchParams.partialComp = 0;
            useAllPatches = 1;
            
            if (exist('shiftVol'))
                firstGuess = shiftVol+1-1;
            else
                firstGuess = single(zeros([4 size(imgVol,2) size(imgVol,3) size(imgVol,4)]));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %carry out the 3D patchMatch
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            
            shiftVol = spatio_temporal_patch_match_mex( imgVol, imgVol,...
                patchMatchParams,firstGuess,occVolDilate,occVolDilate);
         %   if (exist('stop_and_debug.txt'))
         %       keyboard;
         %   end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %carry out the reconstruction
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tRec=tic;
            if (exist('featurePyramid','var') && exist('shiftVol'))
                [imgVol,gradX,gradY,normGradX,normGradY] = reconstruct_video_and_features_mex(imgVol,occVol,...
                    shiftVol,patchMatchParams,...
                    sigmaColour,useAllPatches,reconstructionType);%
                patchMatchParams.gradX = single(gradX); patchMatchParams.gradY = single(gradY); patchMatchParams.normGradX = single(normGradX); patchMatchParams.normGradY = single(normGradY);
            else
                [imgVol] = reconstruct_video_mex(imgVol,occVol,...
                    shiftVol,patchMatchParams,...
                    sigmaColour,useAllPatches,reconstructionType);%
            end
            reconstructionTime=toc(tRec);
            X=sprintf('Reconstruction time %f',reconstructionTime); disp(X)
            iterationNb = iterationNb+1;
        else
            init_coarsest_level;
        end
        pp=pp+1;
        %we have finished the initialisation : make sure that both
        %the patchMatch AND the reconstruction consider that
        %everything is known (no more partial patch comparisons)
        patchMatchParams.partialComparison = 0;
        useAllPatches = 1;  patchMatchParams.useAllPatches = 1;
        
        %calculate the residual to see if we terminate
        residual = sum(abs(imgVolIterMinusOne(:) - imgVol(:)))/(single(3*sum(occVol(:)>0)));
        X=sprintf('Residual : %f',residual); disp(X)
        
    end
    
    
    beep;
    
    if (ii>1)
        imgVol = imgVolFine;%[];%
        %interpolate the shift volume
        shiftVol = single(interpolate_disp_field(shiftVol,imgVol,1/scaleStep, patchSize,'nearest'));
    end
    
    if (ii==1)
        t2 = toc(t1)
        imgVol = reconstruct_video_mex(imgVol,occVol,...
            shiftVol,patchMatchParams,sigmaColour,useAllPatches,1);
        occInds = find(occVol>0);
        energy = sum(shiftVol(occInds + 3*prod(sizeImgVol(1:2))));
        imgOut = imgVol;
        shiftVolOut = shiftVol;
        
        % return;
    end
    %erase some of the video structures
    occVolDilate = single([]);
    occVolPatchMatch = single([]);
    imgVolFine = single([]);
    
    patchMatchParams.partialComp = 0;
    useAllPatches = 1;
end
imgVolOut=imgVol;
%Write to output Image:
cd(rootDataset);
fileOut=strcat(file,num2str(patchSizeX));
if(~exist(fileOut,'dir'))
    mkdir(fileOut);
end;
outputDir= fullfile(rootDataset,fileOut);
cd(outputDir);
write_image_volume(imgVolOut,strcat(fileOut,'inpainted'));
%write time to txt file
cd(currDir);
%create Video file
seqImage_to_aviVideo(outputDir,'png',rootDataset,outputVideoFile);
fileID=fopen('time.txt','a+');
fprintf(fileID,'Filename: %s  Patchsize %d Time: %6.2f \n', file ,patchSizeX, t2);
fclose(fileID);




