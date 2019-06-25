function restoreSpatiallyVariantBlur()
    randn('seed',0);
    dataSize=[512 512];
    
    for idx=1:numel(dataSize)
        ranges{idx}=([1:dataSize(idx)]-floor(dataSize(idx)/2)-1)./dataSize(idx);
    end
    [X,Y]=ndgrid(ranges{:});
    R=sqrt(X.^2+Y.^2);
    mu=0.05;
    sigma=0.025;
    amplitude=@(F) exp(-abs(F-mu).^2./(2*sigma^2));
    
    noiseFt=amplitude(R).*(randn(dataSize)+1i*randn(dataSize))/sqrt(2);
    
    fld=fftshift(ifftn(ifftshift(noiseFt)));
    data=abs(fld).^2;
    data=data./(R/.5).^3>median(data(:));
    
    imagesc(data);
end