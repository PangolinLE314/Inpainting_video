%this script creates a video sequence

function[] = create_video(fileName,ext,compression)

    if (nargin <2)
        ext = '.png';
    end

    currentDir = cd;
    if (~exist(strcat(currentDir,'\','New_video')))
        mkdir 'New_video';
    end

    %get all .tiff files
    currentFiles = dir(strcat(fileName,'*'));
    nbFrames = length(currentFiles);
    
    frameStep = 0.1;

    cd 'New_video';
        if (compression)
            vid = VideoWriter(strcat(fileName,'.avi'));
        else
            vid = VideoWriter(strcat(fileName,'.avi'),'Uncompressed AVI');
        end
        open(vid);
    cd ..

    for ii=1:(frameStep):nbFrames
        currFileName = currentFiles(min(round(ii),nbFrames)).name;
        imgTemp = imread(currFileName);
        %imgTemp = imgTemp(1:380,1:580,:);
        currFrame.cdata = normalise(double(imgTemp));
        if (ndims(currFrame.cdata) == 2 )
            currFrame.cdata = repmat(currFrame.cdata,[1 1 3]);
        end
        currFrame.colormap = [];
        cd 'New_video';
            writeVideo(vid,currFrame);
        cd ..
        %ii
        %writeVideo(vid,currentFrame);
    end
    cd 'New_video';
    close(vid);
    cd ..
    disp('Video generation finished');
    
end
