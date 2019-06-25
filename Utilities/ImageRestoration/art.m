% [x,errors,xNorms] = art(A, b, nbIterations, x0, relaxationParameter)
%
% Algebraic Reconstruction Technique or Kaczmarz iteration
%
function [x, errors, xNorms] = art(A, b, nbIterations, x0, relaxationParameter)
    if nargin<1 || isempty(A),
        nbEqn = 1000;
        nbVar = nbEqn;
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
    
    compareEfficiency = false;
    
    wantPerformanceFigures = nargout~=1 && ~compareEfficiency;
    
    % Compare to MEX C code implementation
    if compareEfficiency,
        % Convert from sparse to full matrix
        A = full(A);
        b = full(b);
        wantPerformanceFigures = false;
        
        startTime = cputime();
        x = ARTReconstructionMex(A, b, relaxationParameter, nbIterations, x0);
        elapsedTimeMEX = cputime()-startTime;
        logMessage('C Calculation time: %0.6fs, error %d', [elapsedTimeMEX norm(b-A*x)]);
    end
    
    % Initialization.
    if wantPerformanceFigures,
        errors = zeros(1, nbIterations);
        xNorms = zeros(1, nbIterations);
    end
    
    if compareEfficiency,
        startTime = cputime();
    end
    
    % precalculate the equation sqd weights
    weights = relaxationParameter./sum(abs(A).^2, 2);
    
    % replace with transpose for efficiency
    A = A.';
    
    % do the Kaczmarz iteration
    x = x0;
    for itIdx = 1:nbIterations,
        for eqnIdx = 1:nbEquations, % should actually be chosen randomly with a distribution proportional to equationSqdNorms
            V = A(:, eqnIdx);
            x = x + (b(eqnIdx)-V'*x)*weights(eqnIdx)*V;
        end
        if wantPerformanceFigures,
            errors(itIdx) = norm(b-A.'*x);
            xNorms(itIdx) = norm(x);
        end
    end
    
    if compareEfficiency,
        elapsedTime = cputime()-startTime;
        logMessage('Matlab Calculation time: %0.6fs, error %d',[elapsedTime norm(b-A.'*x)]);
    end
    clear A;
    
    if nargout == 0 && wantPerformanceFigures,
        close all;
        fig = figure();
        axs(1) = subplot(1,2,1);
        semilogy([1:nbIterations], errors); title('error');
        xlabel('iteration'); ylabel('error');
        axs(2) = subplot(1,2,2);
        semilogy([1:nbIterations], xNorms); title('xNorm');
        xlabel('iteration'); ylabel('xNorm');
        
        linkaxes(axs, 'x');
    end
    if nargout < 1,
        clear x;
    end
end