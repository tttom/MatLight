function output = j2m(input)
%J2M Converts Java object to MATLAB object 
%   Use this function to retrieve structures from java objects.
%   Provided by Mathworks support

    % The first item are the field names
    fields = cellstr(char(input(1)));

    % The second is the data
    data = input(2);
    data = data(1);
    
    % For each field1
    for i=1:length(fields)
        
        field = fields{i};
        
        % Get the data
        d = data(i);
        
        % Based on data assign value directly, perform a cast or even call 
        % this function recursively if there are nested structs 
        
        switch class(d) 
            
            case 'java.lang.Object[]' % Another structure
            
                output.(field) = j2m(d);    
            
            case 'java.lang.String' % A string
                    
                output.(field) = char(d);    
                
            otherwise % Something else which does not need conversion
                
                output.(field) = d;    
            
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Filename:    j2m.m
%
% Author:      Mathworks
%
% Description:
%   This is a MATLAB script that can be used to retrieve object such as 
%   structures from a Java object in the MATLAB environment.
%
%	To call this function:
%		
%   Type 'j2m(input)' in the MATLAB command window where 'input' is the
%   Java object.
%
%   Ensure that the location of this file is in your MATLAB Path.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

