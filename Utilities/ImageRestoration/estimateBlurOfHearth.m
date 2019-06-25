function estimateBlurOfHearth()
    logMessage('Loading 3D image...');
    dataCube=load('C:\Users\tvettenb\Documents\surealism\matlab\ImageRestoration\data\heart_sample5_fc_FEP_5x_new_ex473_em513_s0_r1_o0_rPos1_.mat','dataCubeGood');
    dataCube=dataCube.dataCubeGood;

    dataSize=size(dataCube); dataSize(end+1:3)=1;
    minSubImageSize=[128 128 8];
    
    nbCells=min(dataSize,ceil(2*dataSize./minSubImageSize));
    logMessage('Estimating blur size and SNR...');
    [sigma amplitude noiseLevel xfRange yfRange sigmaFt noiseLevelFt relFitError]=estimateBlurAndSNR(dataCube,[],[],nbCells,minSubImageSize);
    
    SNR=amplitude./noiseLevelFt;
    
    %
    % Output
    %
    save('blurAndSNREstimateHearth_6b.mat');
    
    close all;
    fig=figure('Position',[100 100 800 600]);
    tmp=sigma;
    tmp(SNR<2.0)=NaN;
    tmp=reshape(tmp,[],size(tmp,3));
    plot(repmat(10.0*(([1:size(tmp,2)]-0.5)*minSubImageSize(3)/2),[size(tmp,1) 1]).',squeeze(tmp(:,:)).','LineWidth',3);
    axs=gca();
    labs(1)=xlabel('z [\mum]'); labs(2)=ylabel('\sigma [\mum]');
    set([axs,labs],'FontSize',14,'FontWeight','bold');
    set(axs,'LineWidth',3);
    set(axs,'YLim',[0 get(axs(1),'YLim')*[0;1]]);
    
    fig=figure('Position',[100 100 800 600]);
    for zIdx=1:size(sigma,3),
        SNRSlice=SNR(:,:,zIdx);
        showSlice=SNRSlice>2;
        sigmaSlice=sigma(:,:,zIdx).*showSlice;
        amplitudeSlice=amplitude(:,:,zIdx).*showSlice;
        noiseLevelSlice=noiseLevel(:,:,zIdx).*showSlice;
        
        set(fig,'Name',sprintf('z=%d',zIdx));
        axs(1)=subplot(2,2,1);
        imagesc(sigmaSlice); title('\sigma [pixel]'); colorbar; axis image;
        axs(2)=subplot(2,2,2);
        imagesc(amplitudeSlice); title('signal amplitude'); colorbar;axis image;
        axs(3)=subplot(2,2,3);
        imagesc(noiseLevelSlice*100); title('noise level [%]'); colorbar; axis image;
        axs(4)=subplot(2,2,4);
        imagesc(SNRSlice); title('SNR'); colorbar;axis image;
    %     subplot(1,3,3); imagesc(img); colorbar; axis image;
        set(axs,'CLim',[0 5])
        linkaxes(axs);
        drawnow();
        pause(0.5);
    end
    
end