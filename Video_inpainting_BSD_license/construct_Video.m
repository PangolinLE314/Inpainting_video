% this script construct video data for video inpainting
rootInputDir='D:\PhD\Dataset\Synthesis\chessboard';
inputExtention='jpg';
outputDir= 'D:\PhD\Results\Videos';
nameVideo='chessboard';

%% creat source video
outputName= strcat(nameVideo,'.avi');
seqImage_to_aviVideo(rootInputDir,inputExtention,outputDir,outputName);

%% creat occlusion video
outputNameOcc=strcat(nameVideo,'Occ.avi');
maskDir=fullfile(rootInputDir,'mask');
imageNames = dir(fullfile(rootInputDir,'mask',strcat('*.',inputExtention)));
if(any(size(imageNames,1)))
    seqImage_to_aviVideo(maskDir,inputExtention,outputDir,outputNameOcc);
else if(exist(fullfile(rootInputDir,'gt.txt')))
        imageNames = dir(fullfile(rootInputDir,strcat('*.',inputExtention)));
        load(fullfile(rootInputDir,'gt.txt'))
        occRec=gt;
        occX=occRec(:,1); occY=occRec(:,2); occX2=occRec(:,3); occY2=occRec(:,4);
        imageNames = {imageNames.name}';
        img = imread(fullfile(rootInputDir,imageNames{1}));
        outputVideo = VideoWriter(fullfile(outputDir,outputNameOcc));
        outputVideo.FrameRate = 30;
        open(outputVideo)
        for i=1:size(imageNames,1)
            Xvec=[occX(i) occX2(i) occX2(i) occX(i) occX(i)];
            Yvec=[occY(i) occY(i) occY2(i) occY2(i) occY(i)];
            mask_temp=poly2mask(Xvec,Yvec, size(img,1), size(img,2));
            img(:,:,1)=uint8(mask_temp*255);
            img(:,:,2)=uint8(mask_temp*255);
            img(:,:,3)=uint8(mask_temp*255);
            writeVideo(outputVideo,uint8(img));
        end
        close(outputVideo);
    end
end
