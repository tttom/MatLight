f%
%
%
function displayDataCube(dataCube,aspectRatio)
    if nargin<2 || isempty(aspectRatio),
        aspectRatio=ones(1,ndims(dataCube));
    end
    
    nbSamples=10^4;
    quantiles=[0.01 0.99];
    
    windowSize=[1024 768];
    
    nbColors=1024;
    
    dataSize=size(dataCube);
    % determine the axes sizes so that the horizontal axes are identical
    % and it fits in the window
    physicalSize=dataSize.*aspectRatio;
    spacing=max(min([.05 .05],30./windowSize),60./windowSize);
    scaling(1) = (1-3*spacing(1))/(physicalSize(2)+physicalSize(3));
    scaling(2) = (1-2*spacing(2))/physicalSize(1);
    displaySize=physicalSize*min(scaling);
    realSpacing=([1 1]-displaySize(1:2)-[0 displaySize(3)])./[2 3];
    axPos1=[realSpacing.*[1 2]+[0 displaySize(3)] displaySize(1:2)];
    axPos2=[realSpacing displaySize([1 3])];
    
    % normalize the intensity
    sampleIndexes=([1:nbSamples]-1+rand(1,nbSamples))/nbSamples;
    sampleIndexes=1+floor(sampleIndexes*numel(dataCube));
    sampling=dataCube(sampleIndexes);
    sampling=sort(sampling);
    levels=sampling(max(1,floor(quantiles*end)));
    normalize=@(x) nbColors*max(0,min(1,(x-levels(1))/diff(levels)));
    
    fig=figure('Position',[50 50 windowSize],'Name',sprintf('3D maximum intensity viewer - %dx%dx%d',dataSize),'NumberTitle','off');
    ax(1)=axes('Units','normalized','Position',axPos1);
    im(1)=image(zeros(dataSize(2),dataSize(1),1));
    hold on;
    ptch(1)=patch(1e6*[-1 1 1 -1].',[0 0 0 0].',[0 0 0 0].','Parent',ax(1));
    hold off;
    ax(2)=axes('Units','normalized','Position',axPos2);
    im(2)=image(zeros(dataSize(3),dataSize(1),1));
    hold on;
    ptch(2)=patch(1e6*[-1 1 1 -1].',[0 0 0 0].',[0 0 0 0].','Parent',ax(2));
    hold off;
    set(ax,'YDir','reverse');
    set(ax,'XLim',[1 dataSize(1)]);
    set(ax(1),'YLim',[1 dataSize(2)]);
    set(ax(2),'YLim',[1 dataSize(3)]);
    xlabel('x','Parent',ax(1)); ylabel('y','Parent',ax(1));
    xlabel('x','Parent',ax(2)); ylabel('z','Parent',ax(2));
    colormap(redHotColorMap(nbColors));
    cb=colorbar('Position',[0.9 0.1 0.025 0.8],'YTick',[1:3],'YTickLabel',{'0.0','0.5','1.0'});
    set(ptch,'EdgeColor',[0 0 1],'EdgeAlpha',0.80,'FaceColor',[0 0 1],'FaceAlpha',0.20);
    % store app data
    setappdata(fig,'dataCube',dataCube);
    setappdata(fig,'ax',ax);
    setappdata(fig,'normalize',normalize);
    setappdata(fig,'downPos',floor(dataSize./2)+1);
    setappdata(fig,'upPos',floor(dataSize./2)+1);
    setappdata(fig,'im',im);
    setappdata(fig,'colorMap',colorMap);
    setappdata(fig,'ptch',ptch);
    setappdata(fig,'cb',cb);
    
    updateFig(fig);
end
function updateFig(fig)
    upPos=getappdata(fig,'upPos');
    downPos=getappdata(fig,'downPos');
    im=getappdata(fig,'im');
    ptch=getappdata(fig,'ptch');
    dataCube=getappdata(fig,'dataCube');
    colorMap=getappdata(fig,'colorMap');
    ax=getappdata(fig,'ax');
    normalize=getappdata(fig,'normalize');
    cb=getappdata(fig,'cb');
    
    dataSize=size(dataCube);
    
    downPos=min(dataSize,max(1,round(downPos)));
    upPos=min(dataSize,max(1,round(upPos)));
    ranges=[downPos;upPos];
    ranges=sort(ranges);
    
    xyImg=max(dataCube(:,:,ranges(1,3):ranges(2,3)),[],3).';
    xyImg=normalize(xyImg);
    xzImg=squeeze(max(dataCube(:,ranges(1,2):ranges(2,2),:),[],2)).';
    xzImg=normalize(xzImg);
    
    set(im(1),'CData',xyImg);
    set(ptch(1),'YData',[downPos(2)*[1 1] upPos(2)*[1 1]].');
    set(im(2),'CData',xzImg);
    set(ptch(2),'YData',[downPos(3)*[1 1] upPos(3)*[1 1]].');
    
    drawnow();

    nbColors=size(colormap(ax(1)),1);
    nbLabels=numel(get(cb,'YTickLabel'));
    set(cb,'YTick',max(1,([0:(nbLabels-1)]./(nbLabels-1))*nbColors));
    
    for axIdx=1:2,
        set(ax(axIdx),'ButtonDownFcn',@(src,evt) markPos(axIdx,src));
    end
    set([ptch im],'HitTest','off');
end

function markPos(axIdx,src)
    clickPos=get(src,'CurrentPoint');
    clickPos=clickPos(1,1:2);
    
    fig=get(src,'Parent');
    downPos=getappdata(fig,'downPos');
    if axIdx==1,
        downPos=[clickPos downPos(3)];
    else
        downPos=[clickPos(1) downPos(2) clickPos(2)];
    end
    setappdata(fig,'downPos',downPos);
    
    set(fig,'WindowButtonMotionFcn',@(fig,evt,n) updatePos(axIdx,src,true));
    set(fig,'WindowButtonUpFcn',@(fig,evt,n) updatePos(axIdx,src,false));
end
function updatePos(axIdx,src,cont)
    clickPos=get(src,'CurrentPoint');
    clickPos=clickPos(1,1:2);
    
    fig=get(src,'Parent');
    upPos=getappdata(fig,'upPos');
    if axIdx==1,
        upPos=[clickPos upPos(3)];
    else
        upPos=[clickPos(1) upPos(2) clickPos(2)];
    end
    setappdata(fig,'upPos',upPos);
    
    if ~cont,
        set(fig,'WindowButtonMotionFcn',[]);
    end
    
    updateFig(fig);
end