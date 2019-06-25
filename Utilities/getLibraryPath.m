%
% Returns the full path to the matlab/lib/ folder
%
function basePath=getLibraryPath()
    functionName=mfilename();
    fullPath=mfilename('fullpath');
    basePath=fullPath(1:end-length(functionName));
end