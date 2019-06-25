function xEst = ARTReconstructionScala(A, b, nbIterations, x0, relaxationParameter)
    if nargin<1 || isempty(A),
        nbEqn = 1000;
        nbVar = 1000;
        A = randn(nbEqn, nbVar);
        clear nbEqn nbVar
    end
    [nbEquations, nbVariables] = size(A);
    if nargin<2 || isempty(b),
        b = A*ones(nbVariables, 1);
    end
    if nargin<3 || isempty(nbIterations),
        nbIterations = 100;
    end
    if nargin<4 || isempty(x0),
        x0 = zeros(nbVariables, 1, 'double');
    end
    if nargin<4 || isempty(relaxationParameter),
        relaxationParameter = 1;
    end
    
    classPath = 'C:\Users\Tom\Documents\labmadrid\matlab\Utilities\ImageRestoration\';
    scalaPath = 'C:\Program Files (x86)\scala\lib\scala-library.jar';
    className = 'Tom.ART';
    
    javaaddpath(classPath);
    javaaddpath(scalaPath);

    artObject = getInstance(className);
    
%     tic;
    xEst = artObject.calculate(A, b, nbIterations, x0, relaxationParameter);
%     toc
%     
%     xRef = art(A, b, nbIterations, x0, relaxationParameter);
%     relDiff = norm(xEst - xRef)/norm(xRef)
    
    if nargout<1,
        clear xEst;
    end
    
%     clear artObject
%     javarmpath(classPath);
%     javarmpath(scalaPath);
    
end

function instance = getInstance(className)
    try
       thisClass = java.lang.Class.forName(className);
    catch Exc
       classLoader = com.mathworks.jmi.ClassLoaderManager.getClassLoaderManager;
       thisClass = classLoader.loadClass(className);
    end
    
    instance = thisClass.newInstance();
end