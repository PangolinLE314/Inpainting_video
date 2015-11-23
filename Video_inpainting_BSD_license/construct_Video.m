% this script construct video data for video inpainting
rootInputDir='D:\PhD\Dataset\Background_SubtractionGT\truck\private_truth\color';
inputExtention='bmp';
outputDir= fullfile(cd,'Content');

%% creat source video
outputName= 'truck.avi';
seqImage_to_aviVideo(rootInputDir,inputExtention,outputDir,outputName);

%% creat occlusion video
imageNames = dir(fullfile(rootInputDir,strcat('*.',inputExtention)));
imageNames = {imageNames.name}';
outputNameOcc='rubikOcc.avi';
outputVideo = VideoWriter(fullfile(outputDir,outputNameOcc));
outputVideo.FrameRate = 30;
open(outputVideo)
load(fullfile(rootInputDir,'gt.txt'))
occRec=gt;
occX=occRec(:,1); occY=occRec(:,2); occX2=occRec(:,3); occY2=occRec(:,4);
for i=1:size(imageNames,1)
    img = imread(fullfile(rootInputDir,imageNames{i}));
    Xvec=[occX(i) occX2(i) occX2(i) occX(i) occX(i)];
    Yvec=[occY(i) occY(i) occY2(i) occY2(i) occY(i)];
    mask_temp=poly2mask(Xvec,Yvec, size(img,1), size(img,2));
    img(:,:,1)=uint8(mask_temp*255);
    img(:,:,2)=uint8(mask_temp*255);
    img(:,:,3)=uint8(mask_temp*255);
    writeVideo(outputVideo,img);
end
close(outputVideo);