%This function convert from sequence of image to an AVI file
function [] = seqImage_to_aviVideo(inputDir,fileExtention,outputDir,outputName)
imageNames = dir(fullfile(inputDir,strcat('*.',fileExtention)));
imageNames = {imageNames.name}';
outputVideo = VideoWriter(fullfile(outputDir,outputName));
outputVideo.FrameRate = 30;
open(outputVideo)
gtFile=fullfile(inputDir,'gt.txt');
if(exist(gtFile,'file'))
    load(gtFile)
    occX=gt(:,1); occY=gt(:,2); occX2=gt(:,3); occY2=gt(:,4);
end
for ii = 1:length(imageNames)
    img = imread(fullfile(inputDir,imageNames{ii}));
    if(exist(gtFile,'file'))
        Xvec=[occX(ii) occX2(ii) occX2(ii) occX(ii) occX(ii)];
        Yvec=[occY(ii) occY(ii) occY2(ii) occY2(ii) occY(ii)];
        mask=poly2mask(Xvec,Yvec, size(img,1), size(img,2));
        structEl= strel('arbitrary', ones(3,3));
        mask_temp=imerode(mask,structEl);
        mask_temp= mask-mask_temp;
        img(:,:,1)=img(:,:,1)+uint8(mask*255);
        img(:,:,2)=img(:,:,2)+uint8(mask*255);
        img(:,:,3)=img(:,:,3)+uint8(mask*255);
    end
    writeVideo(outputVideo,img);
end
close(outputVideo);
end