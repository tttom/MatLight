%
% [xOptim,err]=splitBregmanMinimize(x0,l1Function,l2Function,minimizeL2,l1SubstituteConvergenceWeight,nbIterations,stopCriterion,wantStatistics,groundTruth,progressCallBack)
%
% Finds the x that minimizes ||l1Function(x)||_1 for ||l2Function(x)||_2^2<epsilon
% where epsilon converges to zero with increasing iteration.
%
% Inputs:
%     x0: a starting point for the solution. This can be a matrix or anything else that can be
%     passed to the l1Function and l2Function as an argument. If
%     groundTruth is specified it should be the same object and
%     subtractable to enable comparison.
%
%     l1Function: a function that takes a matrix or struct of the same size as x0 as input
%         and returns a shrinkable object such as a column vector. The function or
%         method shrink(obj,l1Weight) can be applied to a shrinkable object,
%         resulting in a shrinkable object of the same size but with values u that
%         minimize l1Weight*||joint(u)||_1+||u-obj||_2^2, where joint(u) is a
%         function that joins 1 or more parts of the input using the l2 norm (isotropic shrinking).
%         Instead of a matrix, an instance of Shrinkable can be used as well.
%
%     l2Function: returns a vector to be minimized in the l2-norm.
%
%     minimizeL2: A function with parameters (l1Target,l2TargetOffset, l1SubstConvergenceWeight[, x0, nbIterations(3)])
%             finds the u that minimizes l1SubstConvergenceWeight^2*||l1Function(u)-l1Target||_2^2 + ||l2Function(u)-l2TargetOffset||_2^2
%             When an empty matrix [] is provided non-linear optimization
%             is performed automatically. The number of iterations can
%             be configured by the argument value nbIterations(3).
%             An optional 4th input argument can be taken as an approximate
%             starting value for the optimization algorithm.
%             An optional 5th input argument specifies the desired number
%             of iterations as given by the user in nbIterations(3), default 1.
%
%     l1SubstituteConvergenceWeight:
%             A scalar indicating the weight of the l1 substitute to make sure that l1-norm converges before the l1Substitute
%             Default: 0.10
%
%     nbIterations: A vector with the number of outer and inner iterations to perform.
%                   A third value can indicate the iterations performed in
%                   the default minimizeL2.
%
%     stopCriterion(outerIt, l1Error, l2Error):
%                   a function that returns true when the outer iteration
%                   is to be stopped before nbIterations(1) is reached. It
%                   can accept the input arguments outerIt, l1Error, and
%                   l2Error, where outerIt is the numer of outer iterations
%                   currently performed, and l1Error is the mean absolute
%                   error of l1Function and l2Error is root-mean-square error
%                   of the l2Function.
%                   When equal to the empty list, the iteration is stopped
%                   when the outer iteration reaches nbIterations(1).
%                   The default is [].
%
%     groundTruth: The ideal solution, when provided this is used to
%     track the evolution of the solution.
%
%     progressCallBack: a call back function executed every iteration, taking the arguments:
%           iterationNb: a vector specifying the current outer to inner iteration numbers
%           totalNbIterations: a vector specifying the total number of iterations that will be performed per loop
%           restored: the current optimization result
%           l1Substitute: the substitute variable for the l1 function
%           l1BregmanConstrainers: the Bregman variable to minimize the
%                   difference between the l1 substitute and the l1 function
%           l2BregmanConstrainer: the Bregman variable to minimize the l2 function
%
% Contact Tom Vettenburg <tom.vettenburg@gmail.com> for more information.
% Leave this not in if reusing this code or parts of it. Thanks.
%
% Outputs:
%     xOptim: The result of size identical to x0
%     err: a structure with several matrixes containing the l1 and l2
%     error, as well as its sum.
%
% Examples:
%     % find the peaks in 1D
%     variableSize=[5 1];
%     l1Weight=10;
%     l1SubstWeight=0.5;
%     x0=zeros(variableSize);
%     groundTruth=zeros(variableSize); groundTruth(2:3)=1;
%     observation=groundTruth+0.01*randn(variableSize);
%     l1Function=@(x) l1Weight*(x); % sparsity
%     l2Function=@(x) x-observation; % denoising
%     minimizeL2=[];
%     nbIterations=[50 10];
%     [xOptim err]=splitBregmanMinimize(x0,l1Function,l2Function,minimizeL2,l1SubstWeight,nbIterations,true,groundTruth);
%     itNumbers=1/nbIterations(2)+[0:1/nbIterations(2):(nbIterations(1)-1/nbIterations(2))];
%     semilogy(itNumbers,cat(2,err.groundTruth(:),err.l1(:),err.l2(:),err.l1Bregman(:),err.l2Bregman(:))); hold on;
%     semilogy(cat(1,err.groundTruth(end,:),err.l1(end,:),err.l2(end,:),err.l1Bregman(end,:),err.l2Bregman(end,:)).','*');
%     title('error'); xlabel('iteration'); legend({'\Delta ground truth','\Delta l1 target','\Delta l2 target','l1Bregman','l2Bregman'});
%     figure(); plot(xOptim);
%
%     % Denoise 3D
%     variableSize=[5 5 5];
%     l1Weight=10;
%     l1SubstWeight=1/10;
%     x0=zeros(variableSize);
%     groundTruth=zeros(variableSize); groundTruth(2:4,2:4,1:3)=1;
%     observation=groundTruth+0.10*randn(variableSize);
%     diff3D=@(x) l1Weight*cat(4,diff(x([end,1:end],:,:,1),1,1),diff(x(:,[end,1:end],:,2),1,2),diff(x(:,:,[end,1:end],3),1,3));
%     l1Function=@(x) Shrinkable(diff3D(repmat(x,[1 1 1 3])),3); % sparsity of gradient
%     l2Function=@(x) x-observation; % denoising
%     minimizeL2=[];
%     splitBregmanMinimize(x0,l1Function,l2Function,minimizeL2,l1SubstWeight,[50 1],true,groundTruth);
% 
%     % Deblur and denoise 3D
%     variableSize=[5 5 5];
%     l1Weight=10;
%     l1SubstWeight=1/10;
%     x0=zeros(variableSize);
%     groundTruth=zeros(variableSize); groundTruth(2:3,2:4,1:3)=1;
%     h=@(x) 0.5*x+0.25*circshift(x,[1 0 0])+0.25*circshift(x,-[1 0 0]);
%     observation=h(groundTruth)+0.01*randn(variableSize);
%     diff3D=@(x) l1Weight*cat(4,diff(x([end,1:end],:,:,1),1,1),diff(x(:,[end,1:end],:,2),1,2),diff(x(:,:,[end,1:end],3),1,3));
%     l1Function=@(x) Shrinkable(diff3D(repmat(x,[1 1 1 3])),3); % sparsity of gradient
%     l2Function=@(x) h(x)-observation; % blur
%     minimizeL2=[];
%     splitBregmanMinimize(x0,l1Function,l2Function,minimizeL2,l1SubstWeight,[50 1],true,groundTruth);
%
%     % Deblur and denoise 3D with user L2 solver
%     variableSize=[5 5 5];
%     l1Weight=100;
%     l1SubstWeight=1/10;
%     x0=zeros(variableSize);
%     groundTruth=zeros(variableSize); groundTruth(2:3,2:4,1:3)=1;
%     h=@(x) 0.5*x+0.25*circshift(x,[1 0 0])+0.25*circshift(x,-[1 0 0]);
%     observation=h(groundTruth)+0.00*randn(variableSize);
%     diff3D=@(x) l1Weight*cat(4,diff(x([end,1:end],:,:,1),1,1),diff(x(:,[end,1:end],:,2),1,2),diff(x(:,:,[end,1:end],3),1,3));
%     l1Function=@(x) Shrinkable(diff3D(repmat(x,[1 1 1 3])),3); % sparsity of gradient
%     l2Function=@(x) h(x)-observation; % blur
%     l1FunctionsOnConj=@(x) diff3D(x.components{1}([1,end:-1:2],[1,end:-1:2],[1,end:-1:2],:));
%     diracInput=zeros(variableSize); diracInput(1,1,1)=1;
%     otf=fftn(h(diracInput));
%     diffMatrix=l1Function(diracInput); clear diracInput;
%     fft3=@(x) fft(fft(fft(x,[],1),[],2),[],3);
%     diffMatrixFftSqd=sum(abs(fft3(diffMatrix.components{1})).^2,4); clear diffMatrix;
%     minimizeL2=@(l1Target,l2TargetOffset,l1wt) ifftn((sum(conj(fft3(l1FunctionsOnConj(l1Target*l1wt))),4)*l1wt + ...
%         conj(otf).*fft3(l2TargetOffset+observation))...
%         ./max(diffMatrixFftSqd*(l1wt^2) + abs(otf).^2,eps(1)),'symmetric');
%     splitBregmanMinimize(x0,l1Function,l2Function,minimizeL2,l1SubstWeight,[100 4],true,groundTruth);
%
function [xOptim,err]=splitBregmanMinimize(x0,l1Function,l2Function,minimizeL2,l1SubstituteConvergenceWeight,nbIterations,stopCriterion,wantStatistics,groundTruth,progressCallBack)
    %
    % Initialize defaults for testing
    %
    if nargin<1 && nargout<1
        close all;
    end    
    if nargin<1
        x0=zeros([1 1]*750);
    end
    % Regularization specific
    if nargin<2 || isempty(l1Function)
        if ~isnumeric(x0)
            error('l1Function should be specified if x0 is not a matrix.');
        end
            
        regularizationWeight=100;
        
        diffX=@(x) diff(x([end,1:end],:),1,1);
        diffY=@(x) diff(x(:,[end,1:end]),1,2);
        l1FunctionsParallel=@(x) regularizationWeight*cat(3,diffX(x(:,:,1)),diffY(x(:,:,2)));
        l1Function=@(x) Shrinkable(l1FunctionsParallel(repmat(x,[1 1 2])),2);
        clear diffX diffY regularizationWeight;
        l1FunctionsOnConj=@(x) l1FunctionsParallel(x.components{1}([1,end:-1:2],[1,end:-1:2],:)); % also used by minimizeL2
        clear l1FunctionsParallel;
        
        % diffMatrixFftSqd also required by the default for minimizeL2
        diracInput=x0; diracInput(1,1)=1;
        diffMatrix=l1Function(diracInput);
        diffMatrixFftSqd=sum(abs(fft2(diffMatrix.components{1})).^2,3); % Also required by the default for minimizeL2
        clear diracInput diffMatrix;
        
        %
        % Define the related input arguments now too
        %
        [observation, groundTruth, otf, noiseLevel] = getExampleObservation(size(x0));
        l2Function=@(x) ifft2(otf.*fft2(x),'symmetric')-observation;
        % groundTruth is now know, so it can be used for diagnostics
        minimizeL2=@(l1Target,l2TargetOffset,l1wt) ifft2((sum(conj(fft2(l1FunctionsOnConj(l1Target*l1wt))),3)*l1wt + ...
                conj(otf).*fft2(l2TargetOffset+observation))...
                ./max(diffMatrixFftSqd*(l1wt^2) + abs(otf).^2,eps(1)),'symmetric');
        clear noiseLevel diffMatrixFftSqd otf observation;
    end
    if ~exist('l2Function','var') || isempty(l2Function)
        error('If the l1Function arguments is defined, so must the l2Function argument.');
    end
    if nargin<5 || isempty(l1SubstituteConvergenceWeight)
        l1SubstituteConvergenceWeight=0.10; % Make sure that l1-norm converges before the l1Substitute
    end
    if nargin<6 || isempty(nbIterations)
        nbIterations=[30 3];
    end
    nbIterations(end+1:3)=1;
    if nargin<7 || isempty(stopCriterion)
        stopCriterion=@(outerIter,l1Error,l2Error) false;
    end
    if nargin<8 && ~exist('wantStatistics', 'var')
        wantStatistics=(nargout>1);
    end
    if nargin<9 && ~exist('groundTruth', 'var')
        groundTruth=[];
    end
    if nargin<10 || isempty(progressCallBack)
        progressCallBack=[];
    end
    
    if ~exist('minimizeL2','var') || isempty(minimizeL2)
        minimizeL2=@(l1Target, l2TargetOffset, l1SubstWeight) ...
            genericMinimizeL2(size(x0), nbIterations(3),...
                @(x) l1SubstWeight*(l1Function(x)-l1Target),... % l1SubstWeight under the l2-norm
                @(x) l2Function(x)-l2TargetOffset);
    end
    
    %
    % Do optimization
    %
    % Initialize l1 substitute and Bregman constrainers
    l1Substitute=l1Function(x0)*0; % Get the size(s) of the L1 vector(s)
    l1BregmanConstrainers=l1Substitute;
    l2BregmanConstrainer=l2Function(x0)*0; % Get the size of the L2 vector
    xOptim=x0;
    clear x0
    
    % Initialize an error struct for analysis
    if ~isempty(groundTruth)
        err=struct();
        err.l1=[];
        err.l2=[];
        err.l1Bregman=[];
        err.l2Bregman=[];
        err.groundTruth=[];
    end
    % Constrain the unconstrained problem
    for outerIter = 1:nbIterations(1)
        % Solve unconstrained problem:
        % min_u sum(abs(l1Function(u))) + sum(abs(l2Function(u)-(0-l2BregmanConstrainer)).^2)
        for innerIter = 1:nbIterations(2)
            % Minimize the l2-norms of l2Function and the distance l1Function-l1Substitute:
            % min_u l1SubstituteConvergenceWeight^2*sum(abs(l1Function(u)-(l1Substitute-l1BregmanConstrainers)).^2) + ...
            %                                     sum(abs(l2Function(u)-(0-l2BregmanConstrainer)).^2);
            allMinL2Args={l1Substitute-l1BregmanConstrainers, 0-l2BregmanConstrainer, l1SubstituteConvergenceWeight, xOptim, nbIterations(3)};
            allMinL2Args={allMinL2Args{1:min(end,nargin(minimizeL2))}};
            xOptim=minimizeL2(allMinL2Args{:});
            clear allMinL2Args;
            
            % The weight of the distance l1Function-l1Substitute is increased first by the l1BregmanConstrainer
            
            % Minimize the l2-distance of d to l1Substitute while keeping the l1-norm of d minimized:
            % min_d sum(abs(d)) + l1SubstituteConvergenceWeight^2*sum(abs(l1FunctionCurrentVectors+l1BregmanConstrainers - d).^2)
            l1FunctionCurrentVectors=l1Function(xOptim);
            l1Substitute=shrink(l1FunctionCurrentVectors+l1BregmanConstrainers,1.0/l1SubstituteConvergenceWeight^2);

            % update Bregman parameter for L1 optimization
            l1BregmanConstrainers=l1FunctionCurrentVectors + (l1BregmanConstrainers-l1Substitute);
            
            if ~isempty(progressCallBack)
                try
                    allInputArgs={[outerIter,innerIter],nbIterations, xOptim, l1Substitute, l1BregmanConstrainers, l2BregmanConstrainer};
                    allInputArgs={allInputArgs{1:min(end,nargin(progressCallBack))}};
                    progressCallBack(allInputArgs{:});
                    clear allInputArgs
                catch Exc
                    fprintf('Progress report callback function progressCallBack failed!\n\nError message details:\n%s\n',Exc.getReport('extended'));
                end
            end
            
            l1Error=mean(abs(l1FunctionCurrentVectors(:)));
            clear l1FunctionCurrentVectors;
            
            % Store errors in error struct if requested
            if wantStatistics
                if ~isempty(groundTruth)
                    if ~isstruct(groundTruth)
                        err.groundTruth(innerIter,outerIter)=sqrt(mean(abs(xOptim(:)-groundTruth(:)).^2));
                    else
                        % In case the ground truth is in the form of a struct
                        fieldNames=sort(fields(groundTruth));
                        groundTruthError=0;
                        for fieldIdx=1:numel(fieldNames)
                            groundTruthField=getfield(groundTruth,fieldNames{fieldIdx});
                            xOptimField=getfield(xOptim,fieldNames{fieldIdx});
                            groundTruthError=groundTruthError+sqrt(mean(abs(xOptimField(:)-groundTruthField(:)).^2));
                        end
                        err.groundTruth(innerIter,outerIter)=groundTruthError;
                        clear fieldNames xOptimField groundTruthField groundTruthError
                    end
                else
                    err.groundTruth(innerIter,outerIter)=Inf;
                end
                err.l1(innerIter,outerIter)=l1Error;
                err.l1Bregman(innerIter,outerIter)=sqrt(mean(abs(l1BregmanConstrainers(:)).^2));
                err.l2Bregman(innerIter,outerIter)=sqrt(mean(abs(l2BregmanConstrainer(:)).^2));
            end
        end
        
        % Add errorL2Target back to the observation to iteratively solve the constrained problem
        l2FunctionCurrentVector = l2Function(xOptim);
        l2BregmanConstrainer = l2FunctionCurrentVector + l2BregmanConstrainer;
        
        l2Error = sqrt(mean(abs(l2FunctionCurrentVector(:)).^2));
        clear l2FunctionCurrentVector;
        
        if wantStatistics
            for innerIter = 1:nbIterations(2)
                err.l2(innerIter,outerIter) = l2Error;
            end
        end
        
        % Check if the stop criterion is reached
        if outerIter<nbIterations(1)
            inputArgs = {outerIter,l1Error,l2Error};
            if stopCriterion(inputArgs{1:nargin(stopCriterion)})
                break;
            end
        end
    end
    clear l1BregmanConstrainers l1Substitute l2BregmanConstrainer;
    
    %
    % Output
    %
    if wantStatistics && nargout<2
        figure();
        nbOuterIt=size(err.l2,2);
        itNumbers=1/nbIterations(2)+[0:1/nbIterations(2):(nbOuterIt-1/nbIterations(2))];
        semilogy(itNumbers,cat(2,err.groundTruth(:),err.l1(:),err.l2(:),err.l1Bregman(:),err.l2Bregman(:))); hold on;
        semilogy(cat(1,cat(1,err.groundTruth(end,:),err.l1(end,:),err.l2(end,:),err.l1Bregman(end,:),err.l2Bregman(end,:)).',NaN*ones(1,5)),'*');
        title('error');
        xlabel('iteration');
        legend({'\Delta ground truth','\Delta l1 target','\Delta l2 target','l1Bregman','l2Bregman'});
        fprintf('Difference with ground truth: %0.8f%%\n',100*err.groundTruth(end));
    end
    if nargout<1 && isnumeric(xOptim)
        figure();
        switch(ndims(xOptim))
            case {1,2}
                if any(size(xOptim)==1)
                    plot(xOptim);
                    xlabel('x'); ylabel('intensity');
                else
                    imagesc(real(xOptim));
                    xlabel('x'); ylabel('y');
                end
            otherwise
                subplot(2,2,1);
                imagesc(xOptim(:,:,1,1).');
                xlabel('x'); ylabel('y'); title('z: front'); axis image;
                colorbar();
                subplot(2,2,2);
                imagesc(xOptim(:,:,1+floor(end/2),1).'); axis image;
                xlabel('x'); ylabel('y'); title('z: center');
                colorbar();
                subplot(2,2,3);
                imagesc(real(squeeze(xOptim(:,1,:,1)).')); axis image;
                xlabel('x'); ylabel('z'); title('y: top');
                colorbar();
                subplot(2,2,4);
                imagesc(real(squeeze(xOptim(:,1+floor(end/2),:,1)).')); axis image;
                xlabel('x'); ylabel('z'); title('y: middle');
                colorbar();
        end
        clear 'xOptim';
    end
end


%
% Auxiliary functions in case not all input arguments are given
%

% Get an example image of a given size (padded 750x750px image)
function [observation, groundTruth, otf, noiseLevel] = getExampleObservation(imgSize)
    noiseLevel=0.01;

    groundTruth=1-getTestImage('usaf');
    imgSize(end+1:3)=1;
    groundTruth(1+imgSize(1),1+imgSize(2),1+imgSize(3))=0;
    groundTruth=groundTruth(1:imgSize(1),1:imgSize(2),1:imgSize(3));
    
    randn('seed',0);
    noise=noiseLevel*randn(size(groundTruth));

    psf=fspecial('gaussian',size(groundTruth),5);
    otf=fft2(ifftshift(psf));
    clear psf;
    % otf is also required by the default for l2Function and minimizeL2
    observation=ifft2(otf.*fft2(groundTruth),'symmetric')+noise;
    clear noise;
end

% Optimizes u for sum_i ||fi(u)||_2^2
%
% Input:
%     varSize: the matrix size of the variable u to optimize
%     nbIterations: number of iterations to perform
%     fi(u): one or more functions to be squared and summed
%
function optimX = genericMinimizeL2(varargin)
    varSize = varargin{1};
    if prod(varSize) > 1000
      logMessage('Number of variables is large (%d), provide the minimizeL2 function to make this efficient.', prod(varSize));
    end
    nbIterations = varargin{2};
    functionsToMinimize = varargin(3:end);
    function result = optimFunc(x)
        result = [];
        for argIdx = 1:numel(functionsToMinimize)
            term = functionsToMinimize{argIdx}(reshape(x, varSize));
            result(end+[1:numel(term)]) = term(:);
        end
    end

    optimX = lsqnonlin(@optimFunc, zeros(prod(varSize), 1), [], [], optimset('MaxIter', nbIterations, 'Display', 'off'));
    
    optimX = reshape(optimX,varSize);
end
