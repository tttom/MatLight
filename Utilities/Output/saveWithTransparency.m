% saveWithTransparency(figHandle,fileName)
%
% Saves a figure as a png image with transparency
%
function saveWithTransparency(figHandle,fileName)
    if (~strcmpi(fileName(end-3:end),'.png'))
        fileName=strcat(fileName,'.png');
    end
    % save the original settings
    oldBackGround = get(figHandle,'Color');
    invertHardCopyStatus = get(figHandle,'InvertHardCopy');
    
    set(figHandle,'InvertHardCopy','off');
    
    drawnow();
    % record the image as-is
    colorImg = copyFigToImage(figHandle);
    if any(oldBackGround~=0)
    	% specify an all-black background and record the image
    	set(figHandle,'Color',[0 0 0]);
        blackImg = copyFigToImage(figHandle);
    else
        blackImg = colorImg;
    end
    if any(oldBackGround~=1)
	    % Specify an all-white background and record the image
	    set(figHandle,'Color',[1 1 1]);
        whiteImg = copyFigToImage(figHandle);
    else
        whiteImg = colorImg;
    end
    
    % Calculate the alpha value from the modulation
    alpha=1.0 - double(max(whiteImg-blackImg,[],3))./256;
    
    % Write the image with alpha channel
    imwrite(colorImg,fileName, 'png', 'BitDepth', 16,'Alpha',alpha);
    
    set(figHandle,'Color',oldBackGround);
    set(figHandle,'InvertHardCopy',invertHardCopyStatus);
end

function img = copyFigToImage(figHandle)
%     print(figHandle,'-dpng','-r300',fileName);
%     colorImg = imread(fileName);
    frm = getFrame(figHandle); % may rescale first time
    frm = getFrame(figHandle);
    img = frm.cdata;
end