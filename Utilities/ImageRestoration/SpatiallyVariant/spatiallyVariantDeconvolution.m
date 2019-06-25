% [restoredImg, nbPsfs, coeffsX, coeffsY, undistortedBlurredImg] = spatiallyVariantDeconvolution(blurredImage,psfsImage,nbPsfs,coeffsX,coeffsY)
%
% Undoes distortion and deconvolves blurredImage with the PSFs taken from the
% psfsImage after undistorting it to a rectangular grid.
%
% CoeffsX/Y are matrices with 2D polynomial coefficients on each row for each color-band in the image.
% The vector nbPsfs contains the number of psfs resp. horizontaly and vertically.
%
% Arguments coeffsX, coeffsY, and nbPsfs are optional, but allow a significant speed up.
% Specify coeffsX=[] to indicate no distortion in the input images, nbPsfs
% is still required to know how the psfsImage has to be segmented
%
% Contact Tom Vettenburg <tom.vettenburg@gmail.com> for more information.
% Leave this not in if reusing this code or parts of it. Thanks.
%
% Usage:
%   blurredImage=imread('usafPetal.png');
%   blurredImage=double(blurredImage(:,:,1))./256;
%   psfsImage=imread('grid50pxPetal.png');
%   psfsImage=double(psfsImage(:,:,1))./256;
%   [restoredImg, coeffsX, coeffsY, nbPsfs]=spatiallyVariantDeconvolution(blurredImage,psfsImage)
%    or directly
%   [restoredImg, nbPsfs, coeffsX, coeffsY]=spatiallyVariantDeconvolution('blurredImageFileName.png','psfsImageFileName.png')
%   restoredImg=spatiallyVariantDeconvolution('blurredImageFileName.png','psfsImageFileName.png',nbPsfs,coeffsX,coeffsY)
%
function [restoredImg, nbPsfs, coeffsX, coeffsY, undistortedBlurredImg] = spatiallyVariantDeconvolution(blurredImage,psfsImage,nbPsfs,coeffsX,coeffsY)
    global errorLevel;
    errorLevel=0;
    
    noiseK=10;
    
    if (nargin<1)
        blurredImage='usafPetal.png';
    end
    if (nargin<2)
        psfsImage='grid50pxPetal.png';
    end
    if (ischar(blurredImage))
        blurredImage=double(imread(blurredImage))./256;
    end
    if (ischar(psfsImage))
        psfsImage=double(imread(psfsImage))./256;
    end
    imgSize=size(blurredImage);
    nbBands=imgSize(3);
        
    %Detect distortion unless specified otherwise
    if (nargin<4)
        tic;
        
        %Order of polynomial approximation for distortion
        order=5;
        %Most of the PSF intensity should be contained in a psfD x psfD square around the peak intensity
        psfD=32;
        %Discart all pixels with a value lower than this fraction of the maximum value
        intensityFractionCutOff=0.001;
        
        if (nargin<3 || isempty(nbPsfs))
            [coeffsX, coeffsY, nbPsfs] = findDistortion(psfsImage,order,psfD,intensityFractionCutOff);
        else
            [coeffsX, coeffsY, nbPsfs] = findDistortion(psfsImage,order,psfD,intensityFractionCutOff,nbPsfs);
        end

        logMessage('Distortion detected in %f seconds.',toc);
    end
        
    %Undo distortion if requested
    if (~isempty(coeffsX))
        tic;
        %Undistort both images
        undistortedBlurredImg=undistortImg(blurredImage,coeffsX,coeffsY);%Extend in all directions
        undistortedPsfsImg=undistortImg(psfsImage,coeffsX,coeffsY,0);%Extrapolate with zero
        logMessage('Image and PSF grid corrected for distortion in %f seconds.',toc);
    else
        %Input images are not distored
        undistortedBlurredImg=blurredImage;
        undistortedPsfsImg=psfsImage;
    end
    
    %
    % Do the image restoration now.
    %
    tic;
    %Scale PSFs to on-average unity
    undistortedPsfsImg=undistortedPsfsImg*prod(nbPsfs)/(sum(undistortedPsfsImg(:))/size(undistortedPsfsImg,3));
    psfSize=[size(undistortedPsfsImg,1) size(undistortedPsfsImg,2)]./nbPsfs([2 1]); %The (fractional) partition size
    deconvSize=2.^ceil(log2(2.*psfSize));
    psfOffset=deconvSize./2-floor(psfSize./2);
    %Extend the image edges
    undistortedBlurredImg=extendImageSmoothly(undistortedBlurredImg,deconvSize);
    
    %
    %Calculate weighting to calculate the mean properly at edges
    %
    weighting=zeros(size(undistortedBlurredImg,1),size(undistortedBlurredImg,2));
    %Calculate the patch weights
    yPart=repmat([1:deconvSize(1)].'-1,[1 deconvSize(2)]);
    xPart=repmat([1:deconvSize(2)]-1,[deconvSize(1) 1]);
    wPart=max(0,psfSize(2)-abs(xPart-(1+deconvSize(2)/2))).*max(0,psfSize(1)-abs(yPart-(1+deconvSize(1)/2)));
    wPart=wPart./max(wPart(:));
    %Loop over all partitions to calculate the weighting
    for (y=1:nbPsfs(2))
        for (x=1:nbPsfs(1))
            weighting(ceil((y-1)*psfSize(1))+[1:deconvSize(1)],ceil((x-1)*psfSize(2))+[1:deconvSize(2)])=weighting(ceil((y-1)*psfSize(1))+[1:deconvSize(1)],ceil((x-1)*psfSize(2))+[1:deconvSize(2)])+wPart;
        end
    end
    %
    %Loop over all partitions to calculate the convolution
    %
    %Calculate the relative spatial frequency for use in the regularization
    [X,Y]=meshgrid([-1:2/deconvSize(2):1-1/deconvSize(2)],[-1:2/deconvSize(1):1-1/deconvSize(1)]);
    fRel=sqrt(X.^2+Y.^2);
    restoredImg=zeros(size(undistortedBlurredImg));
    %Loop over all partitions
    for (y=1:nbPsfs(2))
        for (x=1:nbPsfs(1))
            %Extract the PSF
            psf=undistortedPsfsImg(min(size(undistortedPsfsImg,1),1+ceil((y-1)*psfSize(1)):floor(y*psfSize(1))),min(size(undistortedPsfsImg,2),1+ceil((x-1)*psfSize(2)):floor(x*psfSize(2))),:);
            paddedPsf=zeros([deconvSize size(psf,3)]);
            paddedPsf(psfOffset(1)+[1:size(psf,1)],psfOffset(2)+[1:size(psf,2)],:)=psf;            
            %Extract the part of the image to deconvolve
            imgPart=undistortedBlurredImg(ceil((y-1)*psfSize(1))+[1:deconvSize(1)],ceil((x-1)*psfSize(2))+[1:deconvSize(2)],:);
            %Deconvolve
            otf=fftshift(fft2(ifftshift(paddedPsf)));
            noiseToSignalPowerRatio = calcNoiseToSignalRatioForNaturalSpectra(fRel,otf,noiseK,-1);
            wienerFilter=calcWienerFilter(otf,noiseToSignalPowerRatio,1.0);
            if (size(wienerFilter,3)>1)
                wienerFilter=ifftshift(wienerFilter(:,:,[3 1 2]));%Unshift the first two dimensions
            else
                wienerFilter=ifftshift(wienerFilter);
            end
            if (size(wienerFilter,3)==1 && nbBands>1)
                wienerFilter=repmat(wienerFilter,[1 1 nbBands]);
            end
            restoredImgPart=ifft2(fft2(imgPart).*wienerFilter,'symmetric');
            %Weight and add into the result
            restoredImg(ceil((y-1)*psfSize(1))+[1:deconvSize(1)],ceil((x-1)*psfSize(2))+[1:deconvSize(2)],:)=restoredImg(ceil((y-1)*psfSize(1))+[1:deconvSize(1)],ceil((x-1)*psfSize(2))+[1:deconvSize(2)],:)+restoredImgPart.*repmat(wPart,[1 1 nbBands]);
        end
    end
    restoredImg=restoredImg.*repmat((weighting>0)./max(eps(1),weighting),[1 1 nbBands]); %Calc linear interpolation
    restoredImg=restoredImg([1:imgSize(1)]+deconvSize(1)/2,[1:imgSize(2)]+deconvSize(2)/2,:); %Crop
    restoredImg=restoredImg*mean(blurredImage(:))/mean(restoredImg(:)); %Scale so that the mean intensity stays the same (could be included in prev. steps)
    restoredImg=max(0.0,min(1.0,restoredImg));
    logMessage('Image spatially-variantly deconvolved in %f seconds.',toc);
    
    if (nargout==0)
        showImage(restoredImg);
        clear restoredImg;
    end
end

function extendedImg=extendImageSmoothly(img,sizeIncrement,blurRadius)
    origImgSize=size(img);
    extendedImg=zeros([[size(img,1) size(img,2)]+sizeIncrement , size(img,3)]);
    topLeftOffset=ceil(sizeIncrement./2);
    bottomRightOffset=floor(sizeIncrement./2);
    if (nargin>2 && ~isempty(blurRadius) && blurRadius>0)
        %Zero pad
        extendedImg(topLeftOffset(1)+[1:origImgSize(1)]-1,topLeftOffset(2)+[1:origImgSize(2)]-1,:)=img;
        %Blur with Gaussian
        smoothingKernel=fspecial('gaussian',[1 1].*(2.^ceil(log2(blurRadius*5))),blurRadius);
        smoothingKernel=smoothingKernel./sum(smoothingKernel(:));
        extendedImg=convolve(extendedImg,smoothingKernel);
        %Overwrite the edges by an extension of the 1px blurring of the original edge
        extendedImg=2*extendedImg([topLeftOffset(1)*ones(1,topLeftOffset(1)),topLeftOffset(1)+[1:origImgSize(1)],(end-bottomRightOffset(1)+1)*ones(1,bottomRightOffset(1))],[topLeftOffset(2)*ones(1,topLeftOffset(2)),topLeftOffset(2)+[1:origImgSize(2)],(end-bottomRightOffset(2)+1)*ones(1,bottomRightOffset(2))],:);
        %Overwrite the center bit with the original image
        extendedImg(topLeftOffset(1)+[1:origImgSize(1)]-1,topLeftOffset(2)+[1:origImgSize(2)]-1,:)=img;
    else
        %Extend from a rectanglar edge inside the current image
        croppedEdgeSize=[1 1]*32;
        topLeftOffset=topLeftOffset+croppedEdgeSize;
        bottomRightOffset=bottomRightOffset+croppedEdgeSize;
        origImgSize(1:2)=origImgSize(1:2)-2*croppedEdgeSize;
        %Crop an edge of width croppedEdgeSize
        img=img(1+croppedEdgeSize(1):end-croppedEdgeSize(1),1+croppedEdgeSize(2):end-croppedEdgeSize(2),:);
        %Extend by replication with a width of sizeIncrement./2+croppedEdgeSize
        extendedImg=img([ones(1,topLeftOffset(1)),1:origImgSize(1),end*ones(1,bottomRightOffset(1))],[ones(1,topLeftOffset(2)),1:origImgSize(2),end*ones(1,bottomRightOffset(2))],:);
    end
end

function obj=undistortImg(img,coeffsX,coeffsY,extPolValue)
    %Extend image by one pixel in each direction.
    img=img([1,1:end,end],[1,1:end,end],:);
    if (nargin>=4 && ~isempty(extPolValue))
        %Overwrite extension with extPolValue when specified
        img(:,[1 end],:)=extPolValue;
        img([1 end],:,:)=extPolValue;
    end
    nbBands=size(img,3);
    nbDistortionBands=size(coeffsX,1);
    obj=zeros(size(img));
    for (bandIdx=1:nbDistortionBands)
        %Calculate projection indexes and add one for the image extension.
        I=polyval2(coeffsX(bandIdx,:),[1:size(img,2)],[1:size(img,1)])+1;
        J=polyval2(coeffsY(bandIdx,:),[1:size(img,2)],[1:size(img,1)])+1;
        %Clip to image size
        I=max(1,min(I,size(img,2)-1));
        J=max(1,min(J,size(img,1)-1));
        %Interpolate
        fracI=I-floor(I);
        I=floor(I);
        fracJ=J-floor(J);
        J=floor(J);        
        if (nbDistortionBands~=1 || nbBands==1)
            bandIdxs=bandIdx;
            bandOffset=(bandIdx-1)*size(img,1)*size(img,2);
        else
            bandIdxs=[1:nbBands];
            I=repmat(I,[1 1 nbBands]);
            J=repmat(J,[1 1 nbBands]);
            fracI=repmat(fracI,[1 1 nbBands]);
            fracJ=repmat(fracJ,[1 1 nbBands]);
            bandOffset=repmat(permute([bandIdxs-1],[3 1 2])*size(img,1)*size(img,2),[size(img,1) size(img,2) 1]);
        end
        obj(:,:,bandIdxs)=img(J+(I-1)*size(img,1)+bandOffset).*(1-fracI).*(1-fracJ);
        obj(:,:,bandIdxs)=obj(:,:,bandIdxs)+img(J+(I+1-1)*size(img,1)+bandOffset).*(fracI).*(1-fracJ);
        obj(:,:,bandIdxs)=obj(:,:,bandIdxs)+img(J+1+(I-1)*size(img,1)+bandOffset).*(1-fracI).*(fracJ);
        obj(:,:,bandIdxs)=obj(:,:,bandIdxs)+img(J+1+(I+1-1)*size(img,1)+bandOffset).*(fracI).*(fracJ);
    end
    obj=obj(2:end-1,2:end-1,:);
end

function [coeffsX, coeffsY, gridSize]=findDistortion(psfsImage,order,psfD,intensityFractionCutOff,gridSize)
    global errorLevel;
    
    nbBands=size(psfsImage,3);
    positionsPerBand={};
    minIndex=[0 0];
    maxIndex=[0 0];
    for (bandIdx=1:nbBands)
        singleBandImg=psfsImage(:,:,bandIdx);
        %Find positions of the high intensity peaks
        pos=[];
        [sortedIntensities, pixelIdxs]=sort(singleBandImg(:),1,'descend');
        nbHighIntensityPixels=find(sortedIntensities<intensityFractionCutOff*sortedIntensities(1),1);
        if (~isempty(nbHighIntensityPixels))
            sortedIntensities=sortedIntensities(1:nbHighIntensityPixels-1);
            pixelIdxs=pixelIdxs(1:nbHighIntensityPixels-1);
        end
        mask=ones(size(singleBandImg));
        pos=zeros(length(pixelIdxs),3);
        pos(:,3)=sortedIntensities;
        for (idx=1:length(pixelIdxs))
            [row, col]=ind2sub(size(singleBandImg),pixelIdxs(idx));
            newPos=[row col];
            if (mask(row,col))
                pos(idx,[2 1])=newPos;
                %Blacken out a region around the found position
                mask(max(1,newPos(1)+floor(-psfD/2)):min(size(singleBandImg,1),newPos(1)+floor(psfD/2)),max(1,newPos(2)+floor(-psfD/2)):min(size(singleBandImg,2),newPos(2)+floor(psfD/2)))=0;
            end
        end
        pos=pos(pos(:,1)>0,:);
        if (~isempty(errorLevel) && errorLevel>0)
            %Debug
            hold on;
            scatter(pos(:,1),pos(:,2));
        end

        %Connect peaks
        conns=buildConnections(pos);
        %Filter too short connections and remove those positions altogether
        meanDistSqd=mean(conns(:,3));
        shortConnPos=conns(conns(:,3)<=meanDistSqd/4,1);
        pos(shortConnPos,3)=0;%Set intensity to zero of invalid positions
        conns=conns(conns(:,3)>meanDistSqd/4,:);

        if (~isempty(errorLevel) && errorLevel>0)
            %Debug
            plot([pos(conns(:,1),1) pos(conns(:,2),1)].',[pos(conns(:,1),2) pos(conns(:,2),2)].');
        end

        %Order and number
        nbConns=zeros(1,size(pos,1));
        for (idx=1:size(conns,1))      
            nbConns(conns(idx,1))=nbConns(conns(idx,1))+1;
            nbConns(conns(idx,2))=nbConns(conns(idx,2))+1;
        end
        wellConnected=find(nbConns>=4*2);%At least four bi-directional connections
        [~, centerPosIdx]=min(sum(abs(pos(wellConnected,1:2)-repmat(mean(pos(wellConnected,1:2)),[length(wellConnected) 1])).^2,2));
        centerPosIdx=wellConnected(centerPosIdx);
        imgPos=[0 0]+.5;
        pos(centerPosIdx,4:5)=imgPos;
        pos=addBranches(centerPosIdx,pos,conns,450);
        locatedPosMask=(pos(:,4)~=0);
        if (sum(locatedPosMask)<size(pos,1))
            logMessage('Couldn''t position %d points',size(pos,1)-sum(locatedPosMask));
        end
        %Forget about all the unlocated positions and positions near another position
        pos=pos(locatedPosMask & pos(:,3)>0,:);
        clear conns; %The indexes in conns are incorrect now
        
        if (~isempty(errorLevel) && errorLevel>0)
            %Debug
            scatter(pos(:,1),pos(:,2),'filled');
        end
        
        minIndex=min(minIndex,min(pos(:,4:5))-0.5);
        maxIndex=max(maxIndex,max(pos(:,4:5))-0.5);
        
        positionsPerBand{bandIdx}=pos;
    end
    
    if (nargin<5 || isempty(gridSize))
        %Choose the target grid
        gridOffset=-(minIndex+0.5);
        gridSize=1+maxIndex-minIndex;
    else
        %Use the last band's positions to find the offset
        originCoords=find(pos(:,4)==0.5 & pos(:,5)==0.5);
        if (~isempty(originCoords))
            gridOffset=round(pos(originCoords(1),1:2).*gridSize./[size(psfsImage,2),size(psfsImage,1)] - 0.5) -0.5;
        else
            logMessage('Potential problem: origin coord not marked in grid, assuming it is in the center of the image.');
            gridOffset=floor((gridSize-1)./2)-0.5;
        end
    end
    
    %Fit polynomials to distortion
    coeffsX=zeros(nbBands,(order+1)*(order+2)/2);
    coeffsY=zeros(nbBands,(order+1)*(order+2)/2);
    for (bandIdx=1:nbBands)
        pos=positionsPerBand{bandIdx};
        logMessage('Got %d of %d point-spread functions to interpolate, %d missing.',[size(pos,1) prod(gridSize) prod(gridSize)-size(pos,1)]);
    
        %Calculate the position that the peaks should have been at.
        pos(:,4:5)=pos(:,4:5)+repmat(gridOffset,[size(pos,1) 1]);
        pos(:,4:5)=(pos(:,4:5)+0.5).*repmat([size(psfsImage,2) size(psfsImage,1)]./gridSize,[size(pos,1) 1]);

        %Fit polynomials to the displaced grid
        if (size(pos,1)<(order+1)*(order+2))
            logMessage('Insufficient point-spread functions located.');
        end
        if (isempty(errorLevel) || errorLevel<=0) %Turn off warnings unless the errorLevel is set
            ws = warning('off','all'); 
        end
        
        coeffsX(bandIdx,:) = polyfitweighted2(pos(:,4),pos(:,5),pos(:,1),order, ones(size(pos(:,3))) );
        coeffsY(bandIdx,:) = polyfitweighted2(pos(:,4),pos(:,5),pos(:,2),order, ones(size(pos(:,3))) );
        if (isempty(errorLevel) || errorLevel<=0)
            warning(ws); %Turn warnings back on
        end
    end
end

function pos = addBranches(posIdx,pos,conns,recLim)
    currImgPos=pos(posIdx,4:5);
    neighborConns=conns(conns(:,1)==posIdx,:);
    neighborConns=neighborConns(pos(neighborConns(:,2),4)==0,:);%Remove those that lead to nodes that are already numbered
    angs=neighborConns(:,4);
    nextPosIdxs=[];
    for (connIdx=1:size(neighborConns,1))
        ang=angs(connIdx);
        neighborIdx=neighborConns(connIdx,2);
        if (abs(mod(ang-0+pi,2*pi)-pi)<pi/8)
            pos(neighborIdx,4:5)=currImgPos+[0 -1];
            nextPosIdxs(end+1)=neighborIdx;
        end
        if (abs(mod(ang-pi+pi,2*pi)-pi)<pi/8)
            pos(neighborIdx,4:5)=currImgPos+[0 1];
            nextPosIdxs(end+1)=neighborIdx;
        end
        if (abs(mod(ang-pi/2+pi,2*pi)-pi)<pi/8)
            pos(neighborIdx,4:5)=currImgPos+[-1 0];
            nextPosIdxs(end+1)=neighborIdx;
        end
        if (abs(mod(ang+pi/2+pi,2*pi)-pi)<pi/8)
            pos(neighborIdx,4:5)=currImgPos+[1 0];
            nextPosIdxs(end+1)=neighborIdx;
        end
    end
    if (recLim>0)
        for (nextPosIdxsIdx=1:length(nextPosIdxs))
            pos=addBranches(nextPosIdxs(nextPosIdxsIdx),pos,conns,recLim-1);
        end
    end
end

function conns=buildConnections(pos)
    [~, centerPosIdx]=min(sum(abs(pos(:,1:2)-repmat(mean(pos(:,1:2)),[size(pos,1) 1])).^2,2));
    centerPos=pos(centerPosIdx,1:2);
    pos(:,1)=pos(:,1)-centerPos(1);
    pos(:,2)=pos(:,2)-centerPos(2);
    
    numPos=size(pos,1);
    Y=repmat(pos(:,1),[1 numPos]);
    X=repmat(pos(:,2),[1 numPos]);
    distSqd=(X-X.').^2+(Y-Y.').^2;
    ang=atan2(Y-Y.',X-X.');
        
    conns=[];
    for (idx=1:numPos)
        [sortedDists, sortedDistsIdx]=sort(distSqd(idx,:));
        sortedDistsIdx=sortedDistsIdx(1:min(20,length(sortedDistsIdx)));
        highestIntensity=max(pos(sortedDistsIdx,3));%Find the highest intensity in the neighbourhood
        if (pos(idx,3)>0.5*highestIntensity)
            %Connect all points with high intensity unless there is already one in that direction
            distIdx=2;
            newConns=[];
            while (distIdx<length(sortedDistsIdx))
                if (pos(sortedDistsIdx(distIdx),3)>0.5*highestIntensity)
                    newAng=ang(idx,sortedDistsIdx(distIdx));
                    if (isempty(newConns) || min(abs(mod(newConns(:,4)-newAng+pi,2*pi)-pi))>pi/8)
                        newConns(end+1,1:4)=[idx,sortedDistsIdx(distIdx) sortedDists(distIdx) newAng];
                    end
                end
                distIdx=distIdx+1;
            end
            conns(end+[1:size(newConns,1)],1:4)=newConns;
        end
    end
end










