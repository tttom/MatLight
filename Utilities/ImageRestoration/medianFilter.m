% filteredData=medianFilter(data)
%
% Applies a 3x3 median filter to the 2D input data
%
function data=medianFilter(data)
    dataSize=size(data);
    [shiftsX shiftsY]=ndgrid([-1:1],[-1:1]);
    shifts=cat(2,shiftsX(:),shiftsY(:));
    for zIdx=1:prod(dataSize(3:end)),
        filtered=zeros([dataSize(1:2),size(shifts,1)]);
        for shiftIdx=1:size(shifts,1),
            filtered(:,:,shiftIdx)=circshift(data(:,:,zIdx),shifts(shiftIdx,:));
        end
        data(:,:,zIdx)=median(filtered,3);
    end
    data=reshape(data,dataSize);
end