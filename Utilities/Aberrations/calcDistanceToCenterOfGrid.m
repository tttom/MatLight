% distance=calcDistanceToCenterOfGrid(dataSize,gridPitch)
%
% Determines the distance to the center in a regular plaid grid of
% size dataSize and with a given grid pitch: gridPitch 
%
function distance=calcDistanceToCenterOfGrid(dataSize,gridPitch)
    if nargin<1,
        gridPitch=[.5 1];
    end
    if nargin<2,
        dataSize=[100 150];
    end

    nbDims=numel(gridPitch);
    distance=zeros(dataSize,'single');
    for dimIdx=1:nbDims,
        rng=([0:dataSize(1)-1]-floor(dataSize(1)/2))*gridPitch(dimIdx);
        distance=distance+repmat(rng(:).^2,[1 dataSize(2:end)]);
        distance=shiftdim(distance,1);
        dataSize=dataSize([2:end, 1]);
    end
    distance=sqrt(distance);
    
    if nargout<1,
        close all;
        figure;
        imagesc(distance);
        colorbar();
        clear distance;
    end
end