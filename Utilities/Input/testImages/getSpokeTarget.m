% spokeTarget=getSpokeTarget(imgSize,nbSpokes,radius)
%
% Calculates an antialiased spoke target by oversampling 16 times
%
function spokeTarget=getSpokeTarget(imgSize,nbSpokes,radius)
    if nargin<1 || isempty(imgSize)
        imgSize=512;
    end
    if length(imgSize)<2
        imgSize(2)=imgSize(1);
    end
    if nargin<3
        nbSpokes=36;
    end
    if nargin<4
        radius=1;
    end
    
    spokeTarget=zeros(ceil(imgSize(1)/2)); % Calculate on a quarter of its size only
    %For antialiasing without using too much memory
    oversampling=16;
    for dX=-.5:1/(oversampling-1):.5
        for dY=-.5:1/(oversampling-1):.5
            [X,Y]=meshgrid([-1:2/(imgSize(1)+1):-2/(imgSize(1)+1)]+dX*2/(imgSize(1)+1),[-1:2/(imgSize(2)+1):-2/(imgSize(2)+1)]+dY*2/(imgSize(2)+1));
            pixelAngle=atan2(X,Y); % Measured clockwise from the top
            R2=X.^2+Y.^2;
            spokeTarget=spokeTarget+(R2<=radius^2).*(mod(floor(nbSpokes*pixelAngle./pi+0.5),2)<1);
        end
    end
    spokeTarget=spokeTarget/oversampling^2;
    
    %Expand quadrant to full image
    oddNbSpokesPerHalf=mod(nbSpokes/2,2);
    if oddNbSpokesPerHalf,
        spokeTarget(end+1:end+end,:)=flipud(spokeTarget);
    else
        spokeTarget(end+1:end+end,:)=rot90(spokeTarget,1);
    end
    spokeTarget(:,end+1:end+end)=rot90(spokeTarget,2);

    spokeTarget=spokeTarget(1:imgSize(1),1:imgSize(2)); % Crop for odd sizes be cause we rounded up
end