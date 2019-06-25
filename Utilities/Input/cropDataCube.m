%[dataCube xRange yRange zRange]=cropDataCube(dataCube, xLim,yLim,zLim, xRange,yRange,zRange)
%
function [dataCube xRange yRange zRange]=cropDataCube(dataCube, xLim,yLim,zLim, xRange,yRange,zRange)
    dataCubeSize=size(dataCube);
    if (nargin<4 || isempty(xRange))
        xRange=[1:dataCubeSize(1)];
    end
    if (nargin<5 || isempty(yRange))
        yRange=[1:dataCubeSize(2)];
    end
    if (nargin<6 || isempty(zRange))
        zRange=[1:dataCubeSize(3)];
    end
    
    xRangeSel=xRange>=xLim(1) & xRange<=xLim(end);
    yRangeSel=yRange>=yLim(1) & yRange<=yLim(end);
    zRangeSel=zRange>=zLim(1) & zRange<=zLim(end);
    
    dataCube=dataCube(xRangeSel,yRangeSel,zRangeSel,:);
    xRange=xRange(xRangeSel);
    yRange=yRange(yRangeSel);
    zRange=zRange(zRangeSel);
end