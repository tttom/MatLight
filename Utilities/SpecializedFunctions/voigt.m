% v=voigt(x,gamma,sigma,integratedValue,x0l,x0g)
%
% Calculates the Voigt distribution at points x.
%
% Example:
%    x=[-100:.001:100];
%    tic;
%    v=voigt(x,1,1,10,50,0);
%    toc
%    plot(x,v);
%    sum(v)
%
function v=voigt(x,gamma,sigma,integratedValue,x0l,x0g)
    if (nargin<1)
        x=[-10:.01:10];
    end
    if (nargin<2)
        gamma=1;
    end
    if (nargin<3)
        sigma=1;
    end
    if (nargin<4)
        integratedValue=1;
    end
    if (nargin<5)
        x0l=0;
    end
    if (nargin<6)
        x0g=0;
    end
    
    maxCalculationLength=2^20;
    
    maxGaussianTails=1e-3;
    maxLorentzianTails=1e-3;
    maxLorentzianStdX=tan((1-maxLorentzianTails)*pi/2);
    maxGaussianStdX=erfinv(1-maxGaussianTails)/sqrt(0.5);
    
    dx=diff(x(1:2));
    nbOutputSamples=length(x);
    
    x0=x0l+x0g; % Calculate as if Gaussian is centered at 0 from now on
    
    nbSamplesWithSignificantValuesLorentzian=max(1,ceil(2*gamma*maxLorentzianStdX/dx));
    nbSamplesWithSignificantValuesGaussian=max(1,ceil(2*sigma*maxGaussianStdX/dx));
    calculationLength=min(nbOutputSamples,nbSamplesWithSignificantValuesLorentzian)+nbSamplesWithSignificantValuesGaussian;
    if (calculationLength>maxCalculationLength)
        logMessage('Required calculation length %d is larger than the upper limit %d, will perform the calculation on the limitted range.',[calculationLength maxCalculationLength]);
        calculationLength=maxCalculationLength;
    end
    % Select a range, xCalc, on which to do the calculation
    if (calculationLength<nbOutputSamples)
        [~,centerIdx]=min(abs(x-x0));
        leftCropped=min(max(0,centerIdx-ceil(calculationLength/2)),nbOutputSamples-calculationLength);
        rightCropped=nbOutputSamples-calculationLength-leftCropped;
    else
        leftCropped=0;
        rightCropped=0;
    end
    xCalc=x(leftCropped+1)+dx*[0:(calculationLength-1)];
    
    if (nbSamplesWithSignificantValuesLorentzian>1 && nbSamplesWithSignificantValuesGaussian>1)
        lorentzian=lorentzianDistribution(xCalc,x0,gamma);
        gaussian=ifftshift(gaussianDistribution(xCalc,xCalc(1+floor(end/2)),sigma));
        % Convolve
        v=ifft(fft(lorentzian).*fft(gaussian),'symmetric');
    else
        if (nbSamplesWithSignificantValuesLorentzian>1)
            v=lorentzianDistribution(xCalc,x0,gamma);
        else
            v=gaussianDistribution(xCalc,x0,sigma);
        end
    end
    
    % Pad if needed
    v=v([ones(1,leftCropped), 1:end, end*ones(1,rightCropped)]);
    v([1:leftCropped, (end-(rightCropped-1)):end])=0;
    % Crop if needed
    v=v(1:nbOutputSamples);
    
    % Change amplitude as requested
    v=integratedValue*v;
    
    if (nargout==0)
        figure();
        plot(x,v);
        logMessage('Integral value=%0.3f',sum(v));
        clear v;
    end
end

% Calculates the Gaussian distribution that integrates to 1 over infinity
% as integrated per sample interval
function v=gaussianDistribution(x,x0,sigma)
    if ~isempty(x)
        dx=diff(x(1:2));
    %     v=dx*exp(-0.5*((x-x0)./sigma).^2)/(sigma*sqrt(2*pi));
        xIntervalEdges=[x x(end)+dx]-dx/2-x0;
        expIntegrated=erf(xIntervalEdges*sqrt(0.5)/sigma);
        v=diff(expIntegrated)/2;
    else
        v=[];
    end
end

% Calculates the Lorentzian distribution that integrates to 1 over infinity
% as integrated per sample interval
function v=lorentzianDistribution(x,x0,gamma)
    if ~isempty(x)
        dx=diff(x(1:2));
    %     v=dx*(gamma/pi)./((x-x0).^2+gamma^2);
        xIntervalEdges=[x x(end)+dx]-dx/2-x0;
        lorentzianIntegrated=atan2(xIntervalEdges,gamma);
        v=diff(lorentzianIntegrated)/pi;
    else
        v=[];
    end
end
