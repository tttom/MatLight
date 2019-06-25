classdef Cache
    % Cache class
    %
    % Constructor argument:
    %     cacheFileName: The location on disk for storing the cache
    %     data. When omitted a default location is used.
    %
    % Usage:
    %     c=Cache();
    %     isfield(c,'test') => false
    %     c('test')='Test';
    %     c(struct('a','A','Pi',pi))=randn(100);
    %
    %     isfield(c,'test') => true
    %     c('test') => 'Test'
    %     size(c(struct('a','A','Pi',pi))) => [100 100]
    %
    %     numel(c) => 2
    %
    properties (SetAccess=private)
        cacheFileName;
    end
    properties (Access=private)
        matFile;
        version=0.2;
    end
    methods
        function cache=Cache(cacheFileName)
            if (nargin<1)
%                 functionName=mfilename();
%                 path=mfilename('fullpath');
%                 path=path(1:end-length(functionName));
                path=tempdir;
                cacheFileName=fullfile(path,'DefaultCache.mat');
            end
            cache.cacheFileName=cacheFileName;
            
            if (exist('matfile','file') || exist('matfile','builtin'))
                cache.matFile=matfile(cache.cacheFileName,'Writable',true);
            else
                cache.matFile=[]; % pre-2011 version of Matlab
            end
            if (~exist(cache.cacheFileName,'file'))
                % Create an empty file
                version=cache.version;
                if (~isempty(cache.matFile))
                    cache.matFile.version=version;
                else
                    save(cache.cacheFileName,'version','-v7.3');
                end
            end
        end
        function cache=subsasgn(cache,S,value)
            key=S.subs;
            if (length(key)==1)
                key=key{1};
            else
                logMessage('Error, only one index can be specified per assignment!');
                return;
            end
            keyHash=cache.hash(key);
            if (~isempty(cache.matFile))
                setfield(cache.matFile,keyHash,value);
            else
                if iscell(value)
                    value={value}; % The outer cell {} will be removed by the struct constructor on the next line
                end
                structWithKeyHash=struct(keyHash,value);
                save(cache.cacheFileName,'-struct','structWithKeyHash','-append');
            end
        end
        function value=subsref(cache,S)
            key=S.subs;
            if (length(key)==1)
                key=key{1};
            else
                logMessage('Error, only one index can be specified per assignment!');
                return;
            end
            if (isfield(cache,key))
                keyHash=cache.hash(key);
                if (~isempty(cache.matFile))
                    value=getfield(cache.matFile,keyHash);
                else
                    structWithKeyHashValue=load(cache.cacheFileName,keyHash);
                    value=getfield(structWithKeyHashValue,keyHash);
                end
            else
                value=[];
            end
        end
        function value=isfield(cache,key)
            keyHash=cache.hash(key);
            kys=keys(cache);
            value=ismember(keyHash,kys);
        end
%         function rmfield(cache,key) % Does not work
%             keyHash=cache.hash(key);
%             rmfield(cache.matFile,keyHash);
%         end
        function len=numel(cache)
            len=length(cache.keys);
        end
        function len=length(cache)
            len=numel(cache);
        end
        function sz=size(cache,key)
            if (nargin<2)
                sz=[1 numel(cache)];
            else
                if ~isempty(cache.matFile)
                    keyHash=cache.hash(key);
                    sz=size(cache.matFile,keyHash);
                else
                    structWithKeyHash=load(cache.cacheFileName,keyHash);
                    sz=size(getfield(structWithKeyHash,keyHash));
                end
            end
        end
        function clear(cache)
            version=cache.version;
            save(cache.cacheFileName,'version','-v7.3');
        end
    end
    methods (Access=private)
        function kys=keys(cache)
            if (~isempty(cache.matFile))
                kys=who(cache.matFile);
            else
                kys=who('-file',cache.cacheFileName);
            end
            kysI=cellfun(@(x) ~isempty(regexp(x,'^var_', 'once')),kys);
            kys={kys{kysI}};
        end
    end
    methods (Access=private,Static)
        function value=hash(key)
            value=['var_',DataHash(key)];
        end
    end
end