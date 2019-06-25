% values=hyperGaussian(X,order,sigma,centerOffset)
%
% Returns the values for a smoothed edge aperture, with a maximum of 1 in
% the center and sloping of to 0 at a distance of sigma from centerOffset,
% high values of order indicate a sharper transition
%
% Input:
%     X: position
%     order: order, default 2
%     sigma: standard deviation width, default 1
%     centerOffset: the offset of the center, default 0
%
function values=hyperGaussian(X,order,sigma,centerOffset)
    if (nargin<1)
        % Testing
        X=[-5:.01:5];
    end
    if (nargin<2 || isempty(order))
        order=2;
    end
    if (nargin<4 || isempty(centerOffset))
        Xnorm=X;
    else
        Xnorm=X-centerOffset;
    end
    if (nargin>=3 && ~isempty(sigma))
        Xnorm=Xnorm./sigma;
    end
    
    Xnorm=abs(Xnorm);
    
    if (~isinf(order))
        values=exp(-Xnorm.^order);
    else
        values=Xnorm(1)*0+(Xnorm<1);
    end
    
    if (nargin==0 && nargout==0)
        % Testing
        plot(X,values);
        xlabel('X');
        ylabel('values');
        clear values;
    end
end