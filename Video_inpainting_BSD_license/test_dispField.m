% save('dispField.mat','shiftVol');
% dispField= permute(shiftVol,[2 3 1 4]);
% dispFieldX= squeeze(dispField(:,:,1,:));
% dispFieldY= squeeze(dispField(:,:,2,:));
% dispFieldT= squeeze(dispField(:,:,3,:));
% dispFieldSSD= squeeze(dispField(:,:,4,:));
% figure
% image(dispFieldX(:,:,2));
% figure
% image(dispFieldY(:,:,2));
% figure
% image(dispFieldT(:,:,2));
% figure
% image(dispFieldSSD(:,:,3));
% dbquit;
% make_all_mex_files;
% script_inpaiting_video;

ps=3;
script_inpainting_initvideo;
ps=5;
script_inpainting_initvideo
ps=7;
script_inpainting_initvideo;
ps=9;
script_inpainting_initvideo;
