% unwrappedPhase=unwrapWithGradients(complexField,gradients,referencePos,posToUnwrap,useSparseMatrix)
% 
% Returns the unwrapped phase of a n-dimensional complex array.
%
% Inputs:
%     complexField:
%           An n-d array of complex values with a phase that needs
%           to be unwrapped.
%     gradients: an n+1 dimensional array of phase gradients, the final
%           dimension stores n-dimensional gradients
%     referencePos:
%           Optional coordinate or index to start the unwrapping.
%           Conversion ind2sub or sub2ind is done as required.
%           By default the first found maximum value is used.
%    posToUnwrap:
%           Optional boolean array of points that should be unwrapped.
%           By default all values are unwrapped.
%    useSparseMatrix:
%           Optional boolean indicating if the sparse matrix method should
%           be used or not.
%           Default: the algorithm guesses what is fastest
%
% Outputs:
%     unwrappedPhase:
%           The unwrapped phase in radians.
%           An array of the same size as complexField of real values equal
%           to angle(complexField) + 2*pi*K, where K is an integer array.
%
% Inspired by the function QualityGuidedUnwrap2D.m by B.S. Spottiswoode,
% which in turn is adapted from:
% C. Ghiglia and M. phaseDifference. Pritt, Two-Dimensional Phase Unwrapping:
% Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
%
% Example:
%     unwrapnd(); % unwraps a random set and displays the result
%     unwrappedPhase=unwrapnd(randn([1 1]*100).*exp(1i*pi*rand([1 1]*100)));
%
function unwrappedPhase = unwrapWithGradients(complexField, gradients, referencePos, posToUnwrap, useSparseMatrix)
    debug=false;
    if nargin<1 || isempty(complexField),
        randn('seed',0); rand('seed',0);
        initialSize=[1 1]*8; testSize=[1 1]*128;
        complexField=randn(initialSize).*exp(1i*pi*rand(initialSize));
        complexField(testSize(1),testSize(2))=0;
        complexField=circshift(complexField,-floor(initialSize./2));
        complexField=fftshift(ifft2(complexField));
        [X,Y]=meshgrid([1:testSize(1)]-1-floor(testSize(1)./2),[1:testSize(2)]-1-floor(testSize(2)./2));
        P=atan2(Y,X);
        complexField=abs(complexField).*exp(4i*P);
        clear X Y P initialSize testSize;
    end
    
    % Pre-calculate some matrices:
    dataSize=size(complexField);
    allIndexes=[1:prod(dataSize)].';
    allSubs=ind2subs(dataSize,allIndexes);
    dimCumProds=[1, cumprod(dataSize(1:end-1))]; dimCumProds=[dimCumProds;-dimCumProds]; %always used like this anyway
    phase=angle(complexField);
    magnitude=abs(complexField);
    alreadyUnwrapped=false(dataSize);
    % The result will be the same size as the input
    unwrappedPhase=zeros(dataSize);
    
    if nargin<2 || isempty(gradients),
        gradients=zeros([prod(dataSize),numel(dataSize)]);
        shiftVec=zeros(1,numel(dataSize)); shiftVec(1)=1;
        for dimIdx=1:numel(dataSize)
            gradients(:,dimIdx)=(2*phase-circshift(phase,shiftVec)-circshift(phase,-shiftVec))./2;
            shiftVec=shiftVec([end,1:end-1]);
        end
        clear shiftVec;
        gradients=reshape(gradients,[dataSize,numel(dataSize)]);
    end
    if nargin<3 || isempty(referencePos),
        [~, referenceI]=max(magnitude(:));
        referenceI=referenceI(1);
        referencePos=ind2subs(dataSize,referenceI);
    end
    % convert reference from subs to index
    if numel(referencePos)<numel(dataSize),
        referenceI=referencePos;
        referencePos=ind2subs(dataSize,referenceI);
    else
        referenceCells=num2cell(referencePos);
        referenceI=sub2ind(dataSize,referenceCells{:});
        clear referenceCells;
    end
    % store the reference value in the results
    alreadyUnwrapped(referenceI)=true;
    unwrappedPhase(referenceI)=phase(referenceI);
    % find and set the adjoining coordinates
    adjoinI=-dimCumProds+referenceI;
    adjoinI(1,referencePos<=1)=-1;
    adjoinI(2,referencePos>=dataSize)=-1;
    adjoinI=adjoinI(adjoinI>0);
    
    if nargin<3 || isempty(posToUnwrap),
        posToUnwrap=true(dataSize);
    end
    
    % determine the total variation of the gradient
    %gradVariability=calcGradVariability(phase);
    %cost=gradVariability./posToUnwrap;
    cost=1./(magnitude.*posToUnwrap);
    
    if nargin<4 || isempty(useSparseMatrix),
        useSparseMatrix=false;
    end
    
    % Pre-determine a unique integer integer cost for each location
    [~, sortResult]=sort(cost(:));
    orderNumbers(sortResult,1)=allIndexes; % make an inverse index: the locations from good to bad
    
    if ~useSparseMatrix, % this is fastest for small simple problems
        [sortedOrderNumbers, sortI]=sort(orderNumbers(adjoinI));
        adjoinI=adjoinI(sortI); % so that sortedOrderNumbers==orderNumbers(adjoinI)
    else % this is faster for large complex problems 
        % index indicates the order, value indicates the index in the data 
        sortedOrderNumbers=sparse([],[],[], prod(dataSize),1,prod(dataSize));
        sortedOrderNumbers(orderNumbers(adjoinI))=adjoinI;
    end
    clear adjoin;

    % As long as there are 
    unprocessedAdjoining=true;
    while unprocessedAdjoining,  %Loop until there are no more adjoining pixels
        if ~useSparseMatrix;
            %Pop the first element of the list
            activeI=adjoinI(1);
            sortedOrderNumbers=sortedOrderNumbers(2:end);
            adjoinI=adjoinI(2:end);
        else
            orderOfActiveI=find(sortedOrderNumbers>0,1,'first');
            activeI=sortedOrderNumbers(orderOfActiveI);
            sortedOrderNumbers(orderOfActiveI)=0;
        end
        
        activeSubs=allSubs(activeI,:);
        
        % Select the next coordinates to analyze
        nextIs=dimCumProds+activeI;
        nextIs(2,activeSubs<=1)=-1; % Mark if over the edge
        nextIs(1,activeSubs>=dataSize)=-1; % This order for compatibility with previous code
        
        bestMagnitude=-1;
        for nextIIdx=1:numel(nextIs),
            nextI=nextIs(nextIIdx);
            if nextI>0, % Check if not over edge
                if alreadyUnwrapped(nextI),
                    if magnitude(nextI)>bestMagnitude,
                        bestMagnitude=magnitude(nextI);
                        dimI=floor((nextIIdx-1)/size(nextIs,1))+1;
                        % calculate the magnitude-weighted average of the current and the next partial derivative
                        dPhase=sum(gradients([activeI nextI]+(dimI-1)*prod(dataSize)).*magnitude([activeI nextI]))./sum(magnitude([activeI nextI]));
                        if dimCumProds(nextIIdx)>0, dPhase=-dPhase; end % if the next one is to the front, interpolated backwards
                        referencePhase=unwrappedPhase(nextI)+dPhase; % obtain the referencePos unwrapped phase
                        % Do the actual unwrapping
                        bestPhaseUnwrapped=referencePhase+mod(phase(activeI)-referencePhase+pi,2*pi)-pi;
                    end
                else % not yet unwrapped
                    % Put the eligible, still-wrapped neighbors of this pixels in the adjoin set
                    if posToUnwrap(nextI),
                        if ~useSparseMatrix, % faster to check all
                            if all(nextI~=adjoinI),
                                nextOrderNumber=orderNumbers(nextI);

                                lessI=sortedOrderNumbers<nextOrderNumber;
                                adjoinI=[adjoinI(lessI); nextI; adjoinI(~lessI)];
                                sortedOrderNumbers=[sortedOrderNumbers(lessI); nextOrderNumber; sortedOrderNumbers(~lessI)];
                            end
                        else % it is faster to use the fact the data is sorted
                            sortedOrderNumbers(orderNumbers(nextI))=nextI;
                        end
                    end
                end
            end
        end
        
        % Store the result and mark as unwrapped
        unwrappedPhase(activeI)=bestPhaseUnwrapped;
        alreadyUnwrapped(activeI)=true;
                
        if ~useSparseMatrix,
            unprocessedAdjoining=~isempty(sortedOrderNumbers);
        else
            unprocessedAdjoining=nnz(sortedOrderNumbers)>0;
        end
        
        if (debug),
            imagesc(unwrappedPhase); colorbar(); axis('tight');
            set(gca,'CLim',[-15 15]);
            drawnow();
        end
    end
    
    if nargout<1 && ndims(complexField)==2,
        figure();
        [X,Y]=meshgrid([1:size(complexField,1)],[1:size(complexField,2)]);
        surf(X,Y,unwrappedPhase,mod(unwrappedPhase+pi,2*pi)-pi); colormap(hsv);
        colorbar();
        
        clear unwrappedPhase;
    end
end

%
% gradVar=calcGradVariability(phase)
%
% Returns an indicator of the variation in phase gradient.
%
function gradVar=calcGradVariability(phase)
    dataSize=size(phase);
    nbDims=ndims(phase);
    
    gradVar=zeros(dataSize);
    for dimIdx=1:nbDims,
        if dataSize(1)>1,
            dPhase=zeros(dataSize);
            dPhase(1:end-1,:)=diff(phase); % last row zeros
            
            centralD=dPhase; centralD(2:end-1,:)=centralD(2:end-1,:)+dPhase(1:end-2,:);
            centralD([1 end],:)=2*centralD([1 end],:)-centralD([2 end-1],:); % extrapolate the gradient
            centralD=mod(centralD+pi,2*pi)-pi; % unwrap phase
            % list all immediate neighbors
            shiftedCentralD=zeros([dataSize, 1+2*nbDims]);
            shiftedCentralD(1:prod(dataSize))=centralD;
            shiftVector=zeros(1,nbDims); shiftVector(1)=1;
            for innerDimIdx=1:nbDims,
                shiftedCentralD([1:prod(dataSize)]+prod(dataSize)*(1+2*(innerDimIdx-1)))=...
                            circshift(centralD,-shiftVector);
                shiftedCentralD([1:prod(dataSize)]+prod(dataSize)*(2+2*(innerDimIdx-1)))=...
                            circshift(centralD,shiftVector);
                % cycle through all dimensions
                shiftVector=shiftVector([2:end,1]);
            end
            % mean and variance with the immediate neighbors
            sumDiff=mean(shiftedCentralD,nbDims+1);
            totalVar=mean(abs(shiftedCentralD).^2,nbDims+1);
            gradVar=gradVar+sqrt(totalVar-sumDiff.^2);
        end
        % cycle through all dimensions
        phase=permute(phase,[2:nbDims, 1]);
        gradVar=permute(gradVar,[2:nbDims, 1]);
        dataSize=dataSize([2:nbDims, 1]);
    end
    
    % Penalize the edges and corners
    for dimIdx=1:nbDims,
        gradVar([1 end],:)=gradVar([1 end],:)*2;
        gradVar=permute(gradVar,[2:nbDims, 1]);
    end
end
