%
% A class to represent variables that can be 'shrunk'.
%
classdef Shrinkable
    properties (SetAccess=private)
        components={};
        sphericalDimensions=[];
    end
    methods
        % Input arguments should be a series of matrix-integer pairs
        function mca=Shrinkable(varargin)
            inputs=varargin;
            while ~isempty(inputs),
                if isa(inputs{1},'Shrinkable'),
                    newComponents=inputs{1}.components;
                    mca.components(end+[1:numel(newComponents)])=newComponents;
                    mca.sphericalDimensions(end+[1:numel(newComponents)])=inputs{1}.sphericalDimensions;
                    clear newComponents;
                    inputs=inputs(2:end);
                else
                    % Do some input checking
                    if numel(inputs)<2,
                        error('Constructor arguments should be a list of Shrinkable objects or matrix-integer pairs.');
                    end
                    if ~isscalar(inputs{2}) || inputs{2}<1 || (inputs{2}-round(inputs{2})>eps(inputs{2})),
                        error('Constructor arguments should be a list of Shrinkable objects or matrix-integer pairs, where the integer indicates the number of spherical dimensions and thus strictly positive.');
                    end
                    dataSize=size(inputs{1});
                    possibleSphericalDims=[1 cumprod(dataSize(end:-1:2))];
                    if ~any(possibleSphericalDims==inputs{2}),
                        error('Constructor arguments should be a list of Shrinkable objects or matrix-integer pairs, where the integer indicates the number of spherical dimensions.');
                    end
                    mca.components(end+1)=inputs(1);
                    mca.sphericalDimensions(end+1)=inputs{2};
                    inputs=inputs(3:end);
                end
            end
        end
        % Minimizes x for the sum of all component sums:
        % l1Weight*sum(sqrt(sum(abs(x).^2,spDims))) + sum(abs(x(:)-l2Offsets(:)).^2) for x
        % where spDims is the number of dimensions in to be joint optimized,
        % and l2Offsets are the current component's values
        function mca=shrink(mca,l1Weight)
            if nargin<2 || isempty(l1Weight),
                l1Weight=1;
            end
            for cellIdx=1:numcell(mca),
                mca.components{cellIdx}=shrink(mca.components{cellIdx},l1Weight,mca.sphericalDimensions(cellIdx));
            end
        end
        function mca=plus(mca,B)
            mca=mca.applyToCells(@(x,y) x+y,mca.components,B.components);
        end
        function mca=minus(mca,B)
            mca=mca.applyToCells(@(x,y) x-y,mca.components,B.components);
        end
        function mca=uminus(mca)
            mca=mca.applyToCells(@(x) -x,mca.components);
        end
        function mca=uplus(mca)
        end
        function mca=times(mca,B)
            mca=mca.applyToCells(@(x,y) x.*y,mca.components,B.components);
        end
        function mca=mtimes(mca,B)
            if isa(mca,'Shrinkable')
                if isa(B,'Shrinkable')
                    mca=mca.applyToCells(@(x,y) x*y,mca.components,B.components);
                else
                    mca=mca.applyToCells(@(x) x*B,mca.components);
                end
            else
                mca=mtimes(B,mca);
            end
        end
        function mca=power(mca,B)
            if isa(B,'Shrinkable')
                mca=mca.applyToCells(@(x,y) x.^y,mca.components,B.components);
            else
                mca=mca.applyToCells(@(x) x.^B,mca.components);
            end
        end
        function mca=mpower(mca,B)
            if isa(B,'Shrinkable')
                mca=mca.applyToCells(@(x,y) x^y,mca.components,B.components);
            else
                mca=mca.applyToCells(@(x) x^B,mca.components);
            end
        end
        function n = numArgumentsFromSubscript(obj, s, indexingContext)
          n = 1;
%            if indexingContext == matlab.mixin.util.IndexingContext.Expression
%               n = 1;
%            else
%               n = length(s(1).subs{:});
%            end
        end
        function varargout = subsref(mca, S)
            if strcmp(S(1).type,'()') && strcmp(S(1).subs,':'),
                values=mca.applyToCells(@(x) x(:).',mca.components);
                values=cell2mat(values.components).';
            else
                % It is a property name instead
                switch(S(1).subs),
                    case 'components'
                        values=mca.components;
                    case 'sphericalDimensions'
                        values=mca.sphericalDimensions;
                end
                if numel(S)>1
                  values=subsref(values,S(2:end));
                end
            end
            varargout{1} = values;
        end
        
        function mca=abs(mca)
            mca=mca.applyToCells(@(x) abs(x),mca.components);
        end
        function mca=conj(mca)
            mca=mca.applyToCells(@(x) conj(x),mca.components);
        end
        
        function nb=numcell(mca)
            nb=numel(mca.components);
        end
        function sm=numel(mca)
            mca=mca.applyToCells(@(x) numel(x(:)),mca.components);
            sm=sum(cell2mat(mca.components));
        end
        function sm=sum(mca)
            mca=mca.applyToCells(@(x) sum(x(:)),mca.components);
            sm=sum(cell2mat(mca.components));
        end
        function mn=mean(mca)
            mn=sum(mca)/numel(mca);
        end
    end
    methods (Access=private)
        function mca=applyToCells(mca,operand,varargin)
            mca.components=cellfun(operand,varargin{:},'UniformOutput',false);
        end
    end
end