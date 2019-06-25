% res=ssurf(X,Y,Z,C)
% 
% A replacement for surf that is more friendly to use.
% Data is converted to double and undersampled if required so the system stays responsive.
%
% Use help surf for more information
%
function res=ssurf(varargin)
    switch(nargin),
        case 1
            X=[];
            Y=[];
            Z=double(varargin{1});
            C=Z;
        case 2
            X=double(varargin{1});
            Y=[];
            Z=double(varargin{2});
            C=Z;
        case 3
            X=double(varargin{1});
            Y=double(varargin{2});
            Z=double(varargin{3});
            C=Z;
        otherwise
            X=double(varargin{1});
            Y=double(varargin{2});
            Z=double(varargin{3});
            C=double(varargin{4});
    end
    % Check if 2D data is specified, otherwise try to handle it gracefully
    if (min(size(Z))<=1),
        if (length(Z)>1),
            if ~isempty(X),
                res=plot(X,Z);
            else
                res=plot(Z);
            end
        else
            logMessage('A scalar value instead of a matrix was specified: %f',Z);
            res=Z;
        end
        if (nargout==0),
            clear('res');
        end
        return;
    end
    % 2D data is fed in
    maxPoints=128;
    
    surfSize=size(Z);
    if size(X,2)==surfSize(2),
        xRange=X(1,:);
    elseif size(X,1)==surfSize(2),
        xRange=X(:,1);
    else
        if ~isempty(X),
            logMessage('The number of elements in X doesn''t correspond to the number of columns in Z, using default value for X.');
        end
        xRange=[1:surfSize(2)];
    end
    if size(Y,1)==surfSize(1),
        yRange=Y(:,1);
    elseif size(Y,2)==surfSize(1),
        yRange=Y(1,:);
    else
        if ~isempty(Y),
            logMessage('The number of elements in Y doesn''t correspond to the number of rows in Z, using default value for X.');
        end
        yRange=[1:surfSize(1)];
    end
    [X,Y]=meshgrid(xRange,yRange);
    
    % Reduce the number of samples if too many
    if(any(surfSize>maxPoints)),
        %Only reduce the number of points, never increase it.
        if (surfSize(2)>maxPoints),
            xStep=(X(end)-X(1))/(maxPoints-1);
        else
            xStep=X(1,2)-X(1,1);
        end
        if (surfSize(1)>maxPoints),
            yStep=(Y(end)-Y(1))/(maxPoints-1);
        else
            yStep=Y(2,1)-Y(1,1);
        end        
        [sX,sY]=meshgrid([X(1):xStep:X(end)],[Y(1):yStep:Y(end)]);
        Z = interp2(X,Y,Z,sX,sY,'*linear');
        C = interp2(X,Y,C,sX,sY,'*linear');
        X=sX;
        Y=sY;
    end
    
    res=surf(X,Y,single(Z),C,varargin{5:end});%Convert to single precission to avoid problems with the Axis settings
    if (nargout==0),
        clear('res');
    end
end