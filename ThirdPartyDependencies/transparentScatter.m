% usage example:
% scatterPoints = transparentScatter(randn(5000,1),randn(5000,1),0.1,0.05);
% set(scatterPoints,'FaceColor',[1,0,0]);
function scatterPoints = transparentScatter(x,y,radii,opacities)
    if nargin<3 || isempty(radii)
        radii=5;
    end
    if nargin<4 || isempty(opacities)
        opacities=1.0;
    end
    
    % Replicate singleton arguments
    if numel(x)==1,
        x=repmat(x,size(y));
    end
    if numel(y)==1,
        y=repmat(y,size(x));
    end
    if numel(radii)==1,
        radii=repmat(radii,size(x));
    end
    if numel(opacities)==1,
        opacities=repmat(opacities,size(x));
    end
    
    opacities=min(1,max(0,double(opacities)));

    defaultColors = get(0,'DefaultAxesColorOrder');
    assert(numel(x)==numel(y) && numel(x)==numel(radii) , 'x, y, radii, and opacities should have the same number of elements or be singletons');
    t=[0:pi/10:2*pi];

    scatterPoints=0*x;
    for diskIdx=1:numel(x),
        scatterPoints(diskIdx) = patch(radii(diskIdx)*sin(t)+ x(diskIdx),radii(diskIdx)*cos(t)+ y(diskIdx),defaultColors(1,:),'edgecolor','none');
        %alpha(scatterPoints(diskIdx),opacities);
        set(scatterPoints(diskIdx), 'FaceAlpha', opacities(diskIdx));
    end
end