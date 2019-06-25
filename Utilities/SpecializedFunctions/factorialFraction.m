% fraction=factorialFraction(nominatorArgs,denominatorArgs)
% 
% Calculates prod_i factorial(A_i) ./ prod_j factorial(B_j), without
% forming the actual factorials.
%
% Examples:
%   factorialFraction(100,98) - factorial(100)/factorial(98)
%   factorialFraction([100 8],[24 90]) - factorial(100)*factorial(8)/(factorial(24)*factorial(90))
%   factorialFraction(2000,1998) - (2000*1999)
%  factorialFraction([2000 9992],[1990 10000])
%
function fraction=factorialFraction(nominatorArgs,denominatorArgs)
    % remove common arguments in nominator and denominator
    [nominatorArgs,denominatorArgs]=removeCommon(nominatorArgs,denominatorArgs);
    % Extract a map of the unique factors
    [nomFactorMap,denomFactorMap]=factProdArgsRemoveCommon(nominatorArgs,denominatorArgs);
    % convert to primes
    nomFactors=factorMapToFactors(nomFactorMap);
    denomFactors=factorMapToFactors(denomFactorMap);
    % and do the same for all the primes factors
    [nomFactors,denomFactors]=removeCommon(nomFactors,denomFactors);
    % the remaining primes need to be multiplied or divided with the result
    
    % Make the division carefully, one by one
    nomIdx=1; denomIdx=1;
    fraction=1;
    while nomIdx<=numel(nomFactors) || denomIdx<numel(denomFactors),
        while (fraction<=1 || denomIdx>numel(denomFactors)) && nomIdx<=numel(nomFactors),
            fraction=fraction*nomFactors(nomIdx);
            nomIdx=nomIdx+1;
        end
        while (fraction>=1 || nomIdx>numel(nomFactors)) && denomIdx<=numel(denomFactors)
            fraction=fraction/denomFactors(denomIdx);
            denomIdx=denomIdx+1;
        end
    end
end
% returns sorted lists with common factors cancelled out
function [A,B]=removeCommon(A,B)
    A=sort(A); B=sort(B);
    
    aIdx=1; bIdx=1;
    uniqueA=true(size(A));
    uniqueB=true(size(B));
    while (aIdx<=numel(A) && bIdx<=numel(B) && A(aIdx)<B(bIdx)),
        while(bIdx<=numel(B) && B(bIdx)<A(aIdx)),
            bIdx=bIdx+1;
        end
        if bIdx<=numel(B),
            % B(bIdx)>=A(aIdx)
            % Mark for keeping or deleting
            uniq=B(bIdx)>A(aIdx);
            uniqueA(aIdx)=uniq;
            uniqueB(bIdx)=uniq;
        end
        aIdx=aIdx+1;
    end
    % keep only the unique entries
    A=A(uniqueA);
    B=B(uniqueB);
end

function allFactors=factors(numbers)
    allFactors=[];
    for numIdx=1:numel(numbers),
        f=factor(numbers(numIdx));
        allFactors(end+[1:numel(f)])=f;
    end
end

% Converts a product of factorials to its factors
function factors=factProdArgsToFactors(args)
    factors=[1:sum(args)];
    for argIdx=2:numel(args),
        rng=[1:args(argIdx)];
        factors(sum(args(1:argIdx-1))+rng)=rng;
    end
end

% returns a two row matrix indicating the number of duplicates
% row 1: last integer with this number of duplicates
% row 2: number of duplicates per integer
function factorMap=factProdArgsToFactorMap(args)
    args=sort(args);
    factorMap=cat(1,args,[numel(args):-1:1]);
end

%
% Takes factorial arguments as input and returns filtered factorMaps
%
function [Ar,Br]=factProdArgsRemoveCommon(A,B)
    A=factProdArgsToFactorMap(A);
    B=factProdArgsToFactorMap(B);
    
    % Reduce the first part of the map partially
    aIdx=1; bIdx=1;
    Ar=[]; Br=[];
    while aIdx<=size(A,2) && bIdx<=size(B,2),
        cIdx=min(A(1,aIdx),B(1,bIdx));
        Ar(1,end+1)=cIdx;
        Ar(2,end)=max(0,A(2,aIdx)-B(2,bIdx));
        Br(1,end+1)=cIdx;
        Br(2,end)=max(0,B(2,bIdx)-A(2,aIdx));
        % increment counters as needed
        aIdx=aIdx+(cIdx==A(1,aIdx));
        bIdx=bIdx+(cIdx==B(1,bIdx));
    end
    % Copy the remainder unaltered
    Ar(:,end+[1:(size(A,2)-aIdx+1)])=A(:,aIdx:end);
    Br(:,end+[1:(size(B,2)-bIdx+1)])=B(:,bIdx:end);
end

% Input:
% row 1: last integer with this number of duplicates
% row 2: number of duplicates per integer
function f=factorMapToFactors(m)
    f=[];
    nextIdx=2;
    for mIdx=1:size(m,2),
        if m(2,mIdx)>0,
            newF=factors([nextIdx:m(1,mIdx)]);
            newF=repmat(newF,[m(2,mIdx) 1]);
            f(end+[1:numel(newF)])=newF(:);
        end
        nextIdx=m(1,mIdx)+1;
    end
end

