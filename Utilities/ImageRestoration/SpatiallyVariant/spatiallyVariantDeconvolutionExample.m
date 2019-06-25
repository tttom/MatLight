% Demonstration of te spatially variant deconvolution. 
%
% Contact Tom Vettenburg <tom.vettenburg@gmail.com> for more information.
% Leave this in if reusing this code or parts of it. Thanks.
%
function spatiallyVariantDeconvolutionExample()
    close all;

    %
    % Define a test recording
    %
    s = RandStream('mt19937ar', 'Seed',1);
        
    pixelPitch=[1 1]*1e-6;
    averageSigma=3*min(pixelPitch);
    stdSigma=200*min(pixelPitch);
    sigmaDistributionSigma=50*min(pixelPitch);
    groundTruth=getTestImage('boats');
    imgSize=size(groundTruth);
    
    noiseLevel=0.01;
    
    %
    % Restoration parameters
    %
    nbIterations=[100 4 1];
    stopCriterion=@(outerIt,l1Err,l2Err) outerIt>4 && l2Err<noiseLevel/2;
    
    l1Weight=1/2;
    l1SubstWeight=1;
    sparsenessConditioning=0.0*l1Weight;
    
    %
    % 
    %
    xRange=([1:imgSize(1)]-1-floor(imgSize(1)./2))*pixelPitch(1);
    yRange=([1:imgSize(2)]-1-floor(imgSize(2)./2))*pixelPitch(2);
    [xfRange,yfRange]=calcFrequencyRanges(xRange,yRange);
    [XF,YF]=ndgrid(xfRange,yfRange);
    RF2=XF.^2+YF.^2;
    
    % Generate a random blur variation
    sigmaDistribution=averageSigma+stdSigma*s.randn(imgSize);
    sigmaDistribution=ifft2(fft2(sigmaDistribution).*exp(-RF2./(2*(1/(2*pi*sigmaDistributionSigma)).^2)));
    % Determine in how many cells we can split the image
    nbPsfs=floor(imgSize./(max(sigmaDistribution(:))*5/min(pixelPitch)));
    
    % Split in cells
    cellSize=imgSize./nbPsfs; % floating point
    xCellPosRange=round(1+cellSize(1)*([1:nbPsfs(1)]-0.5));
    yCellPosRange=round(1+cellSize(2)*([1:nbPsfs(2)]-0.5));
    xCellRange=[1:cellSize(1)]-1-floor(cellSize(1)./2);
    yCellRange=[1:cellSize(2)]-1-floor(cellSize(2)./2);
    
    % Determine the blur at the center of each cell
    selectedSigmas=sigmaDistribution(xCellPosRange,yCellPosRange);
    psfs=gaussian2Dparallel(xCellRange,yCellRange,selectedSigmas/min(pixelPitch));
    
    % Apply a spatially variant blur
    h=@(x) spatiallyVariantConvolution(x,psfs);
    recorded=h(groundTruth);
    recorded=real(recorded)+noiseLevel*s.randn(size(recorded)); % add noise
    
    %
    % Deblur and denoise 2D with user L2 solver
    %
    x0=zeros(imgSize);
    diff2D=@(x) l1Weight*cat(3,diff(x([end,1:end],:),1,1),diff(x(:,[end,1:end]),1,2));
    diff2DTranspose=@(x) -l1Weight*cat(3,diff(x([1:end,1],:,1),1,1),diff(x(:,[1:end,1],2),1,2));
    isotropicGradientReg=@(x) Shrinkable(diff2D(x),2); % sparsity of isotropic gradient in 2D
    sparsenessReg=@(x) Shrinkable(sparsenessConditioning*x,1); % sparsity of fluorescence
    l1Function=@(x) Shrinkable(isotropicGradientReg(x),sparsenessReg(x)); % sparsity of gradient and fluorescence
    l2Function=@(x) h(x)-recorded; % blur
    
    % minimize l1wt^2*||l1Function(x)-l1Target||^2 + ||h(x)-recorded-l2TargetOffset||^2
    hTranspose=@(x) spatiallyVariantConvolution(x,psfs,true);
    
    paddedPsfs=psfs;
    paddedPsfs(2*end,2*end,:,:)=0;
    paddedPsfs=circshift(paddedPsfs,[ceil([size(psfs,1) size(psfs,2)]./2) 0 0]);
    otfs=fft2(paddedPsfs);
    impulse=zeros(size(otfs,1),size(otfs,2)); impulse(1)=1;
    impulseResponse = l1Function(impulse);
    PNSR = sum(abs(fft2(impulseResponse.components{1})).^2,3)+abs(fft2(impulseResponse.components{2})).^2;
    filters=conj(otfs)./(abs(otfs).^2+repmat(PNSR,[1 1 nbPsfs]));
    kernels=ifft2(filters);
        
    function xOptim=minimizeL2(l1Target,l2TargetOffset,l1wt,x0,nbKrylov)
        tolKrylov=noiseLevel*1e-6;
        
        shape=@(x) reshape(x,[size(x0) numel(x)/numel(x0)]);
        function result=AFun(x)
            x=shape(x);
            result1=l1wt*diff2D(x);
            result1b=l1wt*sparsenessConditioning*x;
            result2=real(h(x));
            result=cat(3,result1,result1b,result2);
        end
        function result=AtFun(x)
            x=shape(x);
            result1=l1wt*sum(diff2DTranspose(x(:,:,1:2)),3);
            result1b=l1wt*sparsenessConditioning*x(:,:,3);
            result2=real(hTranspose(x(:,:,4)));
            result=result1(:)+result1b(:)+result2(:);
        end
        rhs1=l1Target.components{1};
        rhs1b=l1Target.components{2};
        rhs2=recorded+l2TargetOffset;
        rhs=cat(3,rhs1,rhs1b,rhs2);
        clear rhs1 rhs1b rhs2
        B=2*real(AtFun(rhs));
        function err=errorAb(x)
            err=AFun(x)-rhs;
            err=sum(abs(err(:)).^2);
        end
        function errGrad=gradErrorAb(x)
            errGrad=2*AtFun(AFun(x)-rhs);
        end
        function Hx = errHessian(x)
            x=shape(x);
            Hx=2*AtFun(AFun(x));
        end
        
        xOptim=spatiallyVariantConvolution(recorded+l2TargetOffset,kernels)...
            +l1wt*(sum(diff2DTranspose(rhs(:,:,1:2)),3)+sparsenessConditioning*rhs(:,:,3));
        
        % find x to minimize ||afun(x)-rhs||_2^2
%         [xOptim,flag,relres,iter,resvec]=bicgstab(@errHessian,B(:), tolKrylov, nbKrylov, [],[], x0(:));
%         grdnt=norm(gradErrorAb(xOptim))/norm(gradErrorAb(x0))

       xOptim=shape(xOptim);
    end
    %[restored,err]=
    restored=splitBregmanMinimize(x0,l1Function,l2Function,@minimizeL2,l1SubstWeight,nbIterations,stopCriterion,true,groundTruth,@(iterNb, totalIt, restored, l1Subst, l1Bregman, l2Bregman) displayIterFunction(iterNb, totalIt, restored, l1Subst, l1Bregman, l2Bregman, xRange,yRange,l1Function,l2Function));
    
    tiledPsfs=reshape(permute(psfs,[1 3 2 4]),[size(psfs,1)*size(psfs,3) size(psfs,2)*size(psfs,4)]);
    
    figure;
    axs(1)=subplot(2,2,1);
    imagesc(xRange*1e6,yRange*1e6,sigmaDistribution); axis image; colorbar(); title('sigmaDistribution')
    xlabel('x  (\mum)'); ylabel('y  (\mum)');
    axs(2)=subplot(2,2,2);
    imagesc(xRange*1e6,yRange*1e6,tiledPsfs); axis image; colorbar(); title('psfs')
    xlabel('x  (\mum)'); ylabel('y  (\mum)');
    axs(3)=subplot(2,2,3);
    imagesc(xRange*1e6,yRange*1e6,recorded); axis image; colorbar(); title('recorded')
    xlabel('x  (\mum)'); ylabel('y  (\mum)');
    axs(4)=subplot(2,2,4);
    imagesc(xRange*1e6,yRange*1e6,real(restored)); axis image; colorbar(); title('restored')
    xlabel('x  (\mum)'); ylabel('y  (\mum)');
    linkaxes(axs);
end

% calculates the Gaussian function on a regular rectangular grid for one or
% more standard deviations, so that it integrates to 1 for each 2D x-y slice:
% all(sum(sum(G))==1)
%
% x/yRange: uniformely increasing sample coordinates
% sRange: a matrix with standard deviations
%
% 
function G=gaussian2Dparallel(xRange,yRange,sRange)
    G=zeros(numel(xRange),numel(yRange),numel(sRange));
    for zIdx=1:numel(sRange)
        GX=gaussian(xRange,sRange(zIdx));
        GY=gaussian(yRange,sRange(zIdx));
        G(:,:,zIdx)=GX.'*GY;
    end
    G=reshape(G,[numel(xRange) numel(yRange) size(sRange)]);
end

% calculates the 1D Gaussian function on a regular rectangular grid so that it
% integrates to 1. sum(G)->1 for -infty < x < +infty 
%
% xRange: uniformely increasing sample coordinates
% sigma: scalar standard deviation of the distribution
%
function G=gaussian(xRange,sigma)
    dx=diff(xRange(1:2));
    xRangeExt=[xRange xRange(end)+dx]-dx/2;
    G=0.5*diff(erf(xRangeExt./(sqrt(2)*sigma)));
end

function displayIterFunction(iterationNb, totalNbIterations, restored, l1Substitute, l1BregmanConstrainers, l2BregmanConstrainer, xRange,yRange,l1Function,l2Function)
    nbUpdates=10;

    finalIteration=prod(totalNbIterations(1:2));
    displayIters=round([1:(finalIteration-1)/(nbUpdates-1):finalIteration]);
    currentIter=iterationNb(1)*totalNbIterations(2)+iterationNb(2);
    if any(currentIter == displayIters)
        l1SubstituteComb=sqrt(sum(abs(l1Substitute.components{1}).^2,4));
        l1BregmanConstrainersComb=sqrt(sum(abs(l1BregmanConstrainers.components{1}).^2,4));
        l1Error=l1Function(restored); l1Error=sqrt(sum(abs(l1Error.components{1}).^2,4));
        l2Error=l2Function(restored);

        %close all;
        fig=figure('Position',[50 50 1024 768],'Color',[1 1 1]);
        set(fig,'Position',[50,50,1024,768],'Name',sprintf('iteration %d (outer), %d (inner)',iterationNb));
        axs(1)=subplot(2,3,1);
        imagesc(xRange*1e6,yRange*1e6,real(squeeze(restored(:,:,1))));
        xlabel('x  (\mum)'); ylabel('y  (\mum)'); title('restored'); axis image;
        colorbar();
        axs(2)=subplot(2,3,2);
        imagesc(xRange*1e6,yRange*1e6,real(squeeze(l1Error(:,:,1))));
        xlabel('x  (\mum)'); ylabel('y  (\mum)'); title('l1 error'); axis image;
        colorbar();
        axs(3)=subplot(2,3,3);
        imagesc(xRange*1e6,yRange*1e6,real(squeeze(l2Error(:,:,1))));
        xlabel('x  (\mum)'); ylabel('y  (\mum)'); title('l2 error'); axis image;
        colorbar();
        axs(4)=subplot(2,3,4);
        imagesc(xRange*1e6,yRange*1e6,real(squeeze(l1SubstituteComb(:,:,1))));
        xlabel('x  (\mum)'); ylabel('y  (\mum)'); title('l1Substitute'); axis image;
        colorbar();
        axs(5)=subplot(2,3,5);
        imagesc(xRange*1e6,yRange*1e6,real(squeeze(l1BregmanConstrainersComb(:,:,1))));
        xlabel('x  (\mum)'); ylabel('y  (\mum)'); title('l1BregmanConstrainers'); axis image;
        colorbar();
        axs(6)=subplot(2,3,6);
        imagesc(xRange*1e6,yRange*1e6,real(squeeze(l2BregmanConstrainer(:,:,1))));
        xlabel('x  (\mum)'); ylabel('y  (\mum)'); title('l2BregmanConstrainer'); axis image;
        colorbar();

        linkaxes(axs);

        drawnow();
    end
end