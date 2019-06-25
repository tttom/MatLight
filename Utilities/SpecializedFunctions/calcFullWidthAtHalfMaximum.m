% fullWidthAtHalfMaximum=calcFullWidthAtHalfMaximum(X,Y,method)
%
% Calculates the full width at half the maximum using various methods.
% 
% X must be a vector with the sample coordinates or a matrix with the
%     coordinates for a sample set in each column. If X is a vector than it
%     must have the same number as elements as the first dimension in Y.
%     The same coordinates are assumed for all sample sets. If it is a
%     matrix it must have the same size as Y.
% Y must be a vector with the sample values or a matrix with the values of
%     a sample set in each column.
% method: (default 'BiasedGaussian')
%          'BiasedGaussian': fits a Gaussian and constant offset to the data.
%          'Gaussian': fits an unbiased Gaussian to the data.
%          'Linear': finds the FWHM starting from the maximum value and
%          finding the 50% points by linear interpolation.
%          'LinearBiased': Same as Linear after subtraction of the minimum
%          value.
%
function fullWidthAtHalfMaximum = calcFullWidthAtHalfMaximum(X,Y,method)
    if (nargin<2 || isempty(Y))
        Y=X;
        X=[1:size(Y,1)];
    end
    if (nargin<3)
        method='BiasedGaussian';
    end

    inputSize=size(Y);
    if (any(size(X)<inputSize))
        X=repmat(X(:),[1 inputSize(2:end)]);
    end
    if (prod(inputSize)>max(inputSize))
        fullWidthAtHalfMaximum=zeros([inputSize(2:end) 1]);
        for curveIdx=1:prod(inputSize(2:end))
            fullWidthAtHalfMaximum(curveIdx)=calcSingleFullWidthAtHalfMaximum(X(:,curveIdx),Y(:,curveIdx),method);
        end
    else
        fullWidthAtHalfMaximum=calcSingleFullWidthAtHalfMaximum(X,Y,method);
    end
end

function fullWidthAtHalfMaximum=calcSingleFullWidthAtHalfMaximum(X,Y,method)
    switch (lower(method))
        case {'gaussian','biasedgaussian'}
            fullWidthAtHalfMaximum=calcSingleGaussianFullWidthAtHalfMaximum(X,Y,method);
        case 'linear'
            fullWidthAtHalfMaximum=calcSingleLinearFullWidthAtHalfMaximum(X,Y);
        case 'biasedlinear'
            fullWidthAtHalfMaximum=calcSingleLinearFullWidthAtHalfMaximum(X,Y-min(Y(:)));
        otherwise
            error('Invalid method specified.');
    end
end
function fullWidthAtHalfMaximum=calcSingleLinearFullWidthAtHalfMaximum(X,Y)
    [maxY, centerI]=max(Y);
    if (maxY<=0)
        fullWidthAtHalfMaximum=Inf;
        return;
    end
    
    Y=Y/maxY;
    
    leftI = centerI - (find(Y(centerI-1:-1:1)<=0.5,1,'first') - 1);
    rightI = centerI + (find(Y(centerI+1:end)<=0.5,1,'first') - 1);
    if (isempty(leftI) || isempty(rightI) || leftI==1 || rightI == length(X))
        fullWidthAtHalfMaximum=Inf;
        return;
    end
    
    % Interpolate linearly
    leftX = X(leftI-1) + (X(leftI)-X(leftI-1)) * (0.5-Y(leftI-1))/(Y(leftI)-Y(leftI-1));
    rightX = X(rightI+1) + (X(rightI)-X(rightI+1)) * (0.5-Y(rightI+1))/(Y(rightI)-Y(rightI+1));
    
    fullWidthAtHalfMaximum = rightX - leftX;
end
function fullWidthAtHalfMaximum=calcSingleGaussianFullWidthAtHalfMaximum(X,Y,method)
    % Or use [curve, goodness] = fit(double(X(:)),double(Y(:)),'gauss8');
    Y=Y./max(Y); % Normalize to avoid the search from stopping to early    
    if (strcmpi(method,'BiasedGaussian'))
        offset=min(Y);
    else
        offset=0;
    end
    magnitude=max(Y)-offset;
%     center=mean(X.*(Y-offset))/mean(Y);
    medianFilteredY=Y(:); medianFilteredY=median([medianFilteredY medianFilteredY([2:end 1]) medianFilteredY([end, 1:end-1])].');
    [~, centerI]=max(medianFilteredY);
    center=X(centerI);
    sigma=sqrt(mean(((X-center).^2).*(Y-offset)));
   
    if (method)
        x0=double([offset magnitude center sigma]);
        [x]=fminsearch(@(x) norm(gaussian(x(1),x(2),x(3),x(4),X)-Y),x0,optimset('Display','none','TolX',1e-9,'TolFun',1e6*max(Y(:))));
    else
        x0=double([magnitude center sigma]);
        [x]=fminsearch(@(x) norm(gaussian(0,x(1),x(2),x(3),X)-Y),x0,optimset('Display','none','TolX',1e-9,'TolFun',1e6*max(Y(:))));
    end
    x=num2cell(x);
    if (method)
        [offset magnitude center sigma]=deal(x{:});
    else
        [magnitude center sigma]=deal(x{:});
    end
    sigma=abs(sigma); % only the absolute value is important really
    
    fullWidthAtHalfMaximum=2*sigma*sqrt(-2*log(.5));
    
%     figure; plot(X,Y,'-',X,gaussian(offset,magnitude,center,sigma,X),':');
end

function Y=gaussian(offset,magnitude,center,sigma,X)
    Y=offset+magnitude*exp(-(X-center).^2/(2*sigma^2));
end