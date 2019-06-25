%
% Calculates the rotational average around the z-axis centered at the x-y-plane.
% The pixels must have equal sides in the x and y dimensions.
%
function img=rotAvg(inputData)
    if nargin<1,
       inputData=repmat(permute(getTestImage('boats'),[1 2 3]),[1 1 10]);
    end
        
    dataSize=size(inputData);
    
    thetas=[-pi:.1:pi];

    [X,Y]=meshgrid([1:dataSize(2)]-1-floor(dataSize(2)/2), [1:dataSize(1)]-1-floor(dataSize(1)/2));
    img=zeros([dataSize([2 3])]);
    for tIdx=1:numel(thetas),
        theta=thetas(tIdx);
        XI=X(1,:)*cos(theta);
        YI=X(1,:)*sin(theta);
        for zIdx=1:dataSize(3),
            img(:,zIdx)=img(:,zIdx)+interp2(X,Y,inputData(:,:,zIdx),XI,YI).';
        end
    end

    img=img./numel(thetas);
end