% [wavefront, residue]=findLowerOrderAberrations(complexField,highestOrder,maxPhaseStd)
%
%
function [wavefront, residue]=findLowerOrderAberrations(complexField,highestOrder,maxPhaseStd)
    if nargin<1,
        relAmplitudes=load('C:\Users\tvettenb\Dropbox\aberrationMeasurement633Laser.mat','relAmplitudes');
        complexField=relAmplitudes.relAmplitudes;
        complexField(isinf(complexField))=1;
%         [X,Y]=ndgrid([-8:7]./8,[-8:7]./8);
%         R2=X.^2+Y.^2;
%         complexField=exp(2i*pi*(R2./.4+Y./3+X./2-2.5*X.^2-5*(X.^3+Y.^3))).*exp(-(R2./(2*.5^2)));
%         complexField=complexField+0.05*(randn(size(complexField))+1i*randn(size(complexField)))./sqrt(2);
    end
    amplitudes=abs(complexField);
    phases=angle(complexField);
    if nargin<2 || isempty(highestOrder),
        highestOrder=4;
    end
    if nargin<3 || isempty(maxPhaseStd),
        maxPhaseStd=3;
    end
    
    [wavefront, residue]=findWavefront(phases,1./amplitudes,highestOrder,maxPhaseStd);
    
    if nargout<1,
        close all;
        figure('Position',[100 100 800 600]);
        axs(1)=subplot(2,2,1);
        showImage(complexField); axis tight;
        title('input field');
        axs(2)=subplot(2,2,2);
        imagesc(phases); colorbar();
        title('input phases');
        axs(3)=subplot(2,2,3);
        imagesc(wavefront); colorbar();
        title('output wavefront');
        axs(4)=subplot(2,2,4);
        showImage(abs(complexField).*exp(1i*residue)); axis tight;
        title('output residue');
        linkaxes(axs);
        
        clear wavefront;
    end
end

function [wavefront, residue]=findWavefront(phases,stdDevs,highestOrder,maxPhaseStd)
    dataSize=size(phases);
    nbDims=ndims(phases);
    dataSize(end+1:nbDims)=1;
    
    wavefront=zeros(dataSize);
    if highestOrder>0,
        % Remove higher order terms first
        for dimIdx=1:nbDims,
            diffs=mod(diff(phases,1,dimIdx)+pi,2*pi)-pi;
            nextStdDevs=circshift(stdDevs,-double([1:nbDims]==dimIdx));
            higherOrderStdDevs=sqrt(stdDevs.^2+nextStdDevs.^2);
            higherOrderStdDevs=removeLastElement(higherOrderStdDevs,dimIdx);
            % recurse
            wavefront1D=findWavefront(diffs,higherOrderStdDevs,highestOrder-1,maxPhaseStd/2);
            wavefront=wavefront+undiff(wavefront1D,dimIdx);
        end
    end
    residue=phases-wavefront;
    % Calculate (weighted with inverse of stdDevs) average phase piston,
    % modulus in an interval centered around the most accurate phase value
    [~,minI]=min(stdDevs(:)); minI=minI(1);
    piston=sum((mod(residue(:)-residue(minI)+pi,2*pi)-pi)./stdDevs(:))*sum(stdDevs(:)) + residue(minI);
    pistonStd=sqrt(sum(abs(mod(residue(:)-piston+pi,2*pi)-pi).^2./stdDevs(:))*sum(stdDevs(:)));
    if pistonStd>maxPhaseStd,
        logMessage('Stddev %0.1f rad for highest order %d.',[pistonStd highestOrder]);
        piston=0;
        wavefront=0*wavefront;
        residue=phases;
    end
    
    wavefront=wavefront+piston;
    
    if nargout>1,
        residue=mod(residue-piston+pi,2*pi)-pi;
    end
end

% remove last element in a specified dimension
function M=removeLastElement(M,dimIdx)
    origSize=size(M);
    nbDims=numel(origSize);
    
    permI=[dimIdx,1:dimIdx-1,dimIdx+1:nbDims];
    
    M=permute(M,permI);
    M=M(1:end-1,:);
    newSize=origSize; newSize(dimIdx)=newSize(dimIdx)-1;
    M=reshape(M,newSize(permI));
    M=ipermute(M,permI);
end

% Reverses the diff operation 
function M=undiff(diffM,dimIdx)
    if nargin<2,
        dimIdx=find(size(diffM)>1,1,'first');
        if isempty(dimIdx),
            dimIdx=find(size(diffM)>0,1,'first');
        end
    end

    origSize=size(diffM);
    nbDims=numel(origSize);
    
    permI=[dimIdx,1:dimIdx-1,dimIdx+1:nbDims];
    M=permute(diffM,permI);
    M=cat(1,zeros(1,prod(origSize(permI(2:end)))),cumsum(M,1));
    newSize=origSize; newSize(dimIdx)=newSize(dimIdx)+1;
    
    M=reshape(M,newSize(permI));
    M=ipermute(M,permI);
end
