function res=minimizeL2Tmp(l1Target,l2TargetOffset, l1wt)
    diff3D=@(x) cat(4,diff(x([end,1:end],:,:,1),[],1),diff(x(:,[end,1:end],:,2),[],2),diff(x(:,:,[end,1:end],3),[],3));
    l1Function=@(x) Shrinkable(diff3D(repmat(x,[1 1 1 3])),3); % sparsity of gradient
    l1ConjFunctions=@(x) diff3D(x.components{1}([1,end:-1:2],[1,end:-1:2],[1,end:-1:2],:));
    
    variableSize=[5 5 5];
    diracInput=zeros(variableSize); diracInput(1,1,1)=1;
    h=@(x) 0.5*x+0.25*circshift(x,[1 0 0])+0.25*circshift(x,-[1 0 0]);
    otf=fftn(h(diracInput));
    groundTruth=zeros(variableSize); groundTruth(2:3,2:4,1:3)=1;
    observation=h(groundTruth);
    
    diffMatrix=l1Function(diracInput); clear diracInput;
    
    fft3=@(x) fft(fft(fft(x,[],1),[],2),[],3);
    diffMatrixFftSqd=sum(abs(fft3(diffMatrix.components{1})).^2,4); clear diffMatrix;
    
    fftSolve=@(l1Target,l2TargetOffset,l1wt) ifftn(conj(fft3(sum(l1ConjFunctions(l1Target*(l1wt^2)),4)*(l1wt^2)) + ...
        otf.*conj(fft3(l2TargetOffset+observation)))...
        ./max(diffMatrixFftSqd*(l1wt^4) + abs(otf).^2,eps(1)),'symmetric');
    
    res=fftSolve(l1Target,l2TargetOffset, l1wt);
%     nom=conj(fft3(sum(l1ConjFunctions(l1Target*(l1wt^2)),4)*(l1wt^2)) + otf.*conj(fft3(l2TargetOffset+observation)));
%     den=max(diffMatrixFftSqd*(l1wt^4) + abs(otf).^2,eps(1));
%     res=ifftn(nom./den,'symmetric');
end
