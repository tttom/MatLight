% RGB=spectrumToRGB(wavelengths,wavelengthIntensities,saturateValues)
%
% Returns the red green and blue component to simulate a given spectrum.
% The number of dimensions of the returned matrix is equal to that of the
% wavelengthIntensities argument, and the dimensions are the same except
% for the last dimension which is always three.
%
% wavelengths can be a vector of wavelengths or a function handle returning the intensity for a given wavelength, in the latter case the argument wavelengthIntensities is ignored.
% If wavelengthIntensities is a matrix, its last dimension should equal the number of wavelengths.
%
% saturateValues (optional, default=true): a boolean to indicate if negative RGB values should be de-saturated to be non-negative
% 
function RGB=spectrumToRGB(wavelengths,wavelengthIntensities,saturateValues)
    wavelengthsSpecified=nargin>=1 && ~isempty(wavelengths);
    if ~wavelengthsSpecified,
%         wavelengths=532e-9; %Nd-YAG, frequency-doubled
%         wavelengths=632.8e-9; %He-Ne
        wavelengths=1e-9*[300:800];
    end
    % From John Walker's http://www.fourmilab.ch/documents/specrend/specrend.c
    wavelengthsXYZTable=[380:5:780]*1e-9;
    if (isa(wavelengths,'function_handle'))
        wavelengthIntensities=wavelengths(wavelengthsXYZTable);
        wavelengths=wavelengthsXYZTable;
    else
        if wavelengthsSpecified,
            if (nargin<2 || isempty(wavelengthIntensities))
                wavelengthIntensities=ones(size(wavelengths));
            end
        else
            wavelengthIntensities=permute(eye(numel(wavelengths)),[3 1 2]);
        end
    end
    if (nargin<3 || isempty(saturateValues))
        saturateValues=true;
    end
    
    % a spectrum must be calculated
    if isa(wavelengthIntensities,'function_handle')
        wavelengthIntensities=wavelengthIntensities(wavelengths);
        % diagonalize
        wavelengthIntensities=spdiags(wavelengthIntensities.',0,numel(wavelengths),numel(wavelengths));
    else
        if numel(wavelengthIntensities)==1,
            % Use the same value for every wavelength
            wavelengthIntensities(end+1:numel(wavelengths))=wavelengthIntensities(end);
            % diagonalize
            wavelengthIntensities=spdiags(wavelengthIntensities.',0,numel(wavelengths),numel(wavelengths));
        end
    end
    
    wavelengths=wavelengths(:);
    
    inputSize=size(wavelengthIntensities);
    nbWavelengths=numel(wavelengths);
    if (nbWavelengths==1 && inputSize(end)~=1)
        inputSize(end+1)=1;
    end
    wavelengthIntensities=reshape(wavelengthIntensities,[],nbWavelengths).';
    
    if (nbWavelengths~=inputSize(end))
        logMessage('spectrumToRGB Error: the number of wavelengths should be equal to the last dimension of the second argument.');
        return;
    end
    
%     % Based on   sRGB = spectrumRGB(wavelengths*1e9);

    [lambdaMatch, xFcn, yFcn, zFcn] = colorMatchFcn('CIE_1964');
    wavelengthsXYZTable=lambdaMatch.'*1e-9;
    XYZTable=[xFcn;yFcn;zFcn].';
    
    XYZPerWavelength=interp1(wavelengthsXYZTable,XYZTable,wavelengths,'pchip',0).'; %XYZ in rows
    XYZ=XYZPerWavelength*wavelengthIntensities; % XYZ in rows, data points in columns
    
    RGB=applycform(XYZ.',makecform('xyz2srgb')).';
            
    %Adjust to intensity of input (TODO, probably needs to be log.-scaled)
    RGB=RGB.*repmat(sum(wavelengthIntensities,1),[3 1]);
    if saturateValues,
        approximated=sum(max(RGB,[],1)>1);
        excess=max(0,sum(wavelengthIntensities,1)-1);
        RGB=RGB+repmat(excess./3,[3 1]);
        RGB=min(1,RGB);
        approximated=sum(approximated);
    else
        approximated=0;
    end
    
    RGB=reshape(RGB.',[inputSize(1:end-1) 3]);
    RGB=full(RGB);
    
    %Display
    if (nargout==0 && ndims(RGB)>=3)
        RGB=RGB./max(eps(1),max(RGB(:)));
    
        %Make sure that it can be displayed as a color.
        RGB=max(0,RGB);
        
        figure;
%         axis off;
        image(wavelengths*1e9,1,RGB);
        xlabel('\lambda [nm]','Interpreter','Tex');
        
        if approximated>0,
            logMessage('Approximated %0.0f%% of the colors!',100*approximated);
        end
        
        clear RGB;
    end
end