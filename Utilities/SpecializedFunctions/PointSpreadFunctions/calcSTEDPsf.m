function calcStedPsf()
    close all;

    % General constants
    NAvogadro=6.0221417930e23;
    
    % Objectives and wavelengths
    excitation={};
    excitation.wavelength=491e-9; % rsGFP 
    excitation.objective.numericalAperture=1.49;
    excitation.objective.immersionRefractiveIndex=1.5;
    excitation.objective.magnification=60;
    excitation.objective.tubeLength=200e-3; %for Mitutoyo the rear focal length is 200mm
    excitation.averagePower=5e-3; % W
    excitation.pulseTime=200e-12; % s
    excitation.repetitionRate=100e6; % Hz
    excitation.efficiency=0.5;
    depletion=excitation;
    depletion.wavelength=650e-9;
    depletion.averagePower=1e-3; % W
    
    detection=excitation;
    detection.wavelength=510e-9; % rsEGFP
      
    pixelPitch=[[1 1]*20e-9 500e-9];
    gridSize=[[1 1]*128 1];
    xRange=([1:gridSize(1)]-floor(gridSize(1)/2)-1)*pixelPitch(1); % light-sheet propagation axis
    yRange=([1:gridSize(2)]-floor(gridSize(2)/2)-1)*pixelPitch(2); % up
    zRange=([1:gridSize(3)]-floor(gridSize(3)/2)-1)*pixelPitch(3); % detection axis
    
    % Sample medium specification
    sample={};
    sample.fluorophore={};
    sample.fluorophore.maxDensity=0.01*1e3*NAvogadro; % m^-3, 0.01*mol/L density 
    sample.fluorophore.extinctionCoefficient=47000*100/NAvogadro; % /fluorophore /m, rsEGFP
    sample.fluorophore.quantumYield=0.36; % rsEGFP
    sample.refractiveIndex=1.33;
    sample.fluorophore.fluorescenceLifeTime=3e-9; % seconds
    % The following is still given per W/m^2 instead of photon/m^2 
    sample.fluorophore.excitationRate=(2e3/0.01^2)/0.02e-3; % for <0.02 ms at 2kW/cm^2 405nm, rsEGFP
    sample.fluorophore.deactivationRate=(0.6e3/0.01^2)/1e-3; % for 1 ms at 0.6kW/cm^2 488nm, rsEGFP
    
    detector.quantumEfficiency=0.3;
    
    % Calculate the spatial irradiance distributions
    excitationPupil=@(U,V) (sqrt(U.^2+V.^2)<=1);
    depletionPupil=@(U,V) (sqrt(U.^2+V.^2)<=1).*exp(1i*atan2(V,U));
%     detectionPupil=@(U,V) (sqrt(U.^2+V.^2)<=1);
    % Calculate the number of photons/pulse/pixel area
    excitationPsf=calcLinearPsf(xRange,yRange,zRange,excitation,excitationPupil,sample.refractiveIndex);
    depletionPsf=calcLinearPsf(xRange,yRange,zRange,depletion,depletionPupil,sample.refractiveIndex);
%     detectionPsf=calcLinearPsf(xRange,yRange,zRange,detection,detectionPupil,sample.refractiveIndex);
%     % Create photon noise for excitation and depletion
%     excitationPsf=poissrnd(excitationPsf);
%     depletionPsf=poissrnd(depletionPsf);

    %Calculate the fluorophore state in the sample
    pixelVolume=prod(pixelPitch);
    % Number of fluorophores in a specific state (dimension 4)
    fluoDist=cat(4,pixelVolume*sample.fluorophore.maxDensity*ones(gridSize),zeros(gridSize),zeros(gridSize));
    excitedElectronsPerPhotonPerMolecule=pixelPitch(3)*sample.fluorophore.extinctionCoefficient*sample.fluorophore.quantumYield;
    relaxedElectronsPerPhotonPerMolecule=excitedElectronsPerPhotonPerMolecule;
    excitationRate=sample.fluorophore.maxDensity*excitedElectronsPerPhotonPerMolecule*excitationPsf;
    spontaneousEmissionRate=1/sample.fluorophore.fluorescenceLifeTime;
    stimulatedEmissionRate=sample.fluorophore.maxDensity*relaxedElectronsPerPhotonPerMolecule*depletionPsf;
    fluoDist=calcStateEvolution(fluoDist,excitationRate,spontaneousEmissionRate,stimulatedEmissionRate,pulseTime);
    stedPsf=fluoDist(:,:,:,2); % those left in the excited state
    clear fluoDist;
    
    % Simulate the detection efficiency
    objectiveSolidAngle=2*pi*(1-sqrt(1-(detection.objective.numericalAperture/detection.objective.immersionRefractiveIndex)^2));
    stedPsf=stedPsf*objectiveSolidAngle/(4*pi)*detection.efficiency*detector.quantumEfficiency;
    
%     % Create photon noise for detection
%     stedPsf=poissrnd(stedPsf);
    
    %% Output
    figure;
    subplot(2,2,1,'FontWeight','bold','FontName','Times','FontSize',18,'LineWidth',3);
    axs(1)=showImage(abs(excitationPsf),-1,xRange*1e6,yRange*1e6); axis equal; title('excitation');
    subplot(2,2,2,'FontWeight','bold','FontName','Times','FontSize',18,'LineWidth',3);
    axs(2)=showImage(abs(depletionPsf),-1,xRange*1e6,yRange*1e6); axis equal; title('depletion');
    subplot(2,2,3,'FontWeight','bold','FontName','Times','FontSize',18,'LineWidth',3);
    axs(3)=showImage(abs(stedPsf),-1,xRange*1e6,yRange*1e6); axis equal; title('STED');
    subplot(2,2,4,'FontWeight','bold','FontName','Times','FontSize',18,'LineWidth',3);
    stedOtf=fftshift(fft2(ifftshift(stedPsf)));
    [XOtf YOtf]=calcOtfGridFromSampleFrequencies(1./pixelPitch(1:2),gridSize(1:2));
    XOtfRange=XOtf(1,floor(end/2)+1:end)/1e6;
    plot(XOtfRange([1 end]),[0 0],'Color',[1 1 1]*.25,'LineWidth',1);
    hold on;
    plot(XOtfRange,[abs(stedOtf(floor(end/2)+1,floor(end/2)+1:end));abs(stedOtf(floor(end/2)+1:end,floor(end/2)+1)).'].','LineWidth',3);
    xlim(XOtfRange([1 end])); ylim([0 max(abs(stedOtf(:)))]);
    title('OTF');
    xlabel('\nu_X  [cycles/\mum]');
    linkaxes(axs);
    
    peakPhotons=max(stedPsf(:))
    totalPhotons=sum(stedPsf(:))
end

% For a unity pupil it returns the probability of finding a photon pases each pixel area
function psf=calcLinearPsf(xRange,yRange,zRange,lightPath,pupilFunctor,refractiveIndex)
    hPlanck=6.6260695729e-34; % m^2*kg/s
    cLight=299792458;

    psf=calcVectorialPsf(xRange,yRange,zRange,lightPath.wavelength,@(X,Y) pupilFunctor(X,Y)/sqrt(2),@(X,Y) 1i*pupilFunctor(X,Y)/sqrt(2),...
        lightPath.objective.numericalAperture,refractiveIndex,lightPath.objective.magnification,lightPath.objective.tubeLength);
    % sum(psf(:)) should be 1
    laserPeakPower=lightPath.averagePower/(lightPath.pulseLength*lightPath.repetitionRate);
    psf=psf*laserPeakPower*lightPath.efficiency; % W
    
    pixelArea=diff(xRange(1:2))*diff(yRange(1:2));
    peakIrradiance=(max(psf(:))/pixelArea)/1e3*0.01^2  % irradiance in kW/cm^2
    
%     psf=psf./max(psf(:));
    photonEnergy=hPlanck*cLight/lightPath.wavelength;
    psf=psf/photonEnergy; % In photons per pixel area
end

function convolvedImg=convolve(img,psf)
    imgFft=fft2(img);
    otf=fft2(ifftshift(psf));
    convolvedImg=ifft2(imgFft.*otf);
end

function fluoDist=calcStateEvolution(fluoDist,excitationRate,spontaneousEmissionRate,stimulatedEmissionRate,timeInterval)
    inputSize=size(fluoDist);
    % excitation pumps from ground to excited only
    % depletion pumps from excited to ground only
%     states={'ground','excited','bleached'};
    nbStates=inputSize(4);
    emissionRate=spontaneousEmissionRate+stimulatedEmissionRate;
    transitionMatrix=zeros([inputSize(1:3) [1 1]*nbStates]);
    transitionMatrix(:,:,:, 2,1)=emissionRate;
    transitionMatrix(:,:,:, 1,2)=excitationRate;
%     transitionMatrix=[0              excitationRate  0;...
%                       emissionRate   0               0;...
%                       0              0               0];
    %Make sure that the sum of the rows is 0:
    for stateIdx=1:nbStates,
        transitionMatrix(:,:,:, stateIdx,stateIdx)=transitionMatrix(:,:,:, stateIdx,stateIdx)-sum(transitionMatrix(:,:,:, stateIdx,:),5);
    end
%     transitionMatrix=transitionMatrix-diag(sum(transitionMatrix,2));
    
    fluoDist=reshape(fluoDist,[],nbStates);
    transitionMatrix=reshape(transitionMatrix,[],nbStates,nbStates);
    
    % dS=S*transitionMatrix * dt, where S is state row vector
    if (timeInterval<Inf)
        % Transient state
        transitionMatrix=transitionMatrix*timeInterval;
        for (voxelIdx=1:size(fluoDist,1))
            fluoDist(voxelIdx,:)=fluoDist(voxelIdx,:)*expm(squeeze(transitionMatrix(voxelIdx,:,:)));
        end
    else
        % Steady state
        % S*transitionMatrix==0
        for (voxelIdx=1:size(fluoDist,1))
            V=null(squeeze(transitionMatrix(voxelIdx,:,:)).');
            V=V(:,1).';
            V=V./sum(V);
            fluoDist(voxelIdx,:)=sum(fluoDist(voxelIdx,:),2)*V;
        end
    end
    
    fluoDist=reshape(fluoDist,inputSize);
end