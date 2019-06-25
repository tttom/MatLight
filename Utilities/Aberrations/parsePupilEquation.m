% [pupilEquationFunctor argumentsUsed]=parsePupilEquation(pupilEquation,X,Y,time,lambda)
%
% Parses a text string that describes a pupil function into an executable
% matlab function and exectutes it if X and Y are specified.
%
% input:
%    pupilEquation: A text string containing the case insensitive variable
%        names: "X", "Y", "Rho", "Phi", "time", and "Lambda", or shorthands
%        "R", "P", "t", "L" for the latter four. The t must be lower case for
%        now to avoid problems with backwards incompatibilities.
%        Only point-wise operations are allowed, hence "*" is converted to
%        ".*", "/" to "./", and "^" to ".^".
%    X and Y: optional matrixes of the same dimensions that define the carthesian
%        coordinates. When specified, the parsed expression is executed
%        at these coordinates and the resulting matrix is returned.
%    time and lambda: scalar values that are only required if the
%        parsed expression contains these variables. The time variables is
%        suggested to be in seconds from an initial time, and the lambda
%        variables is suggested to be the wavelength in meters.
%
% output:
%    If only one input argument is specified, a function handle is returned
%    with four input arguments: X,Y,time,lambda. These arguments are substituted
%    in the parsed expression, and Rho is defined as sqrt(X.^2+Y.^2), while Phi
%    is defined as atan2(Y,X). The arguments X and Y must be matrices of
%    the same dimensions, and time and lambda are scalars.
%    When three or more arguments are supplied to this function, then the
%    function handle is executed for these arguments and a matrix of the
%    size of X and Y is returned instead.
%
%    The second output argument 'argumentsUsed' is a vector of booleans
%    that indicates of respectivelly X, Y, time, and lambda are used.
%
% Contact Tom Vettenburg <tom.vettenburg@gmail.com> for more information.
% Leave this message in if reusing this code or parts of it. Thanks.
%
function [pupilEquationFunctor, argumentsUsed]=parsePupilEquation(pupilEquation,X,Y,time,lambda)
    pupilEquation=strcat('(',strtrim(pupilEquation),')*1.0');
%     pupilEquation(pupilEquation==unicode2native('?'))='R';
%     pupilEquation(pupilEquation==unicode2native('?'))='P';
    pupilEquation=regexprep(pupilEquation,'(^|[^\w])R(HO)?([^\w]|$)','$1sqrt(X.^2+Y.^2)$3','ignorecase');
    pupilEquation=regexprep(pupilEquation,'(^|[^\w])time([^\w]|$)','$1time$2','ignorecase');
    pupilEquation=regexprep(pupilEquation,'(^|[^\w])Th(eta)?([^\w]|$)','$1P$3','ignorecase');
    pupilEquation=regexprep(pupilEquation,'(^|[^\w])(Ph|F)(i)?([^\w]|$)','$1P$4','ignorecase');
    pupilEquation=regexprep(pupilEquation,'(^|[^\w])P([^\w]|$)','$1atan2(Y,X)$2');
    pupilEquation=regexprep(pupilEquation,'(^|[^\w])t([^\w]|$)','$1timeInSeconds$2');
    pupilEquation=regexprep(pupilEquation,'(^|[^\w])L(AMBDA)?([^\w]|$)','$1lambda$3','ignorecase');
    pupilEquation=regexprep(pupilEquation,'([^\.])([\*/\^])','$1.$2','ignorecase');%Force all matrix operations to be element-wise.
    pupilEquation=regexprep(pupilEquation,'\&\&','&'); %Convert logical operators to bit operators
    pupilEquation=regexprep(pupilEquation,'\|\|','|'); %Convert logical operators to bit operators
    
    % Convert U and V to X and Y, respectively
    pupilEquation=regexprep(pupilEquation,'(^|[^\w])U([^\w]|$)','$1X$2');
    pupilEquation=regexprep(pupilEquation,'(^|[^\w])V([^\w]|$)','$1Y$2');
    
    % Check which arguments were used
    argumentsUsed(1)=~isempty(regexp(pupilEquation,'[^\w\d]X[^\w\d]', 'once'));
    argumentsUsed(2)=~isempty(regexp(pupilEquation,'[^\w\d]Y[^\w\d]', 'once'));
    argumentsUsed(3)=~isempty(regexp(pupilEquation,'[^\w\d]timeInSeconds[^\w\d]', 'once'));
    argumentsUsed(4)=~isempty(regexp(pupilEquation,'[^\w\d]lambda[^\w\d]', 'once'));
    % Check if no X, Y would be used
    if (~any(argumentsUsed(1:2)))
        % if so make sure that no scalar is returned
        pupilEquation=strcat(pupilEquation,'*ones(size(X))');
    end
    
    %Do some quick tests before we break the program: 
    try
        pupilEquationFunctor = @(X,Y,timeInSeconds,lambda) eval(pupilEquation);
        testResult = pupilEquationFunctor([-.1 .2; -.1 .2],[-.1 -.1; .2 .2],0);
    catch Exc
        logMessage('Couldn''t parse pupil phase equation ''%s'', blocking the pupil.',pupilEquation);
        pupilEquationFunctor = @(X,Y,time,lambda) zeros(size(X));
        argumentsUsed = [false false false false];
    end
    
    % Execute the function if X and Y are specified and return a matrix instead.
    if (nargin>=3)
        if (nargin>=4 && isempty(time))
            time=0;
        end
        if (nargin>=5 && isempty(lambda))
            lambda=1e-6;
        end
        switch nargin
            case 1+2
                pupilEquationFunctor=pupilEquationFunctor(X,Y);
            case 1+3
                pupilEquationFunctor=pupilEquationFunctor(X,Y,time);
            otherwise
                pupilEquationFunctor=pupilEquationFunctor(X,Y,time,lambda);
        end
    end
end