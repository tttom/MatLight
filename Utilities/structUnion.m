% C=structUnion(A,B)
%
% Combines arrays or structs A and B recursively where leafs (strings and
% scalars) of B have priority. This function is useful in combination with
% the loadjson function.
%
function C=structUnion(A,B)
    if (isstruct(A))
        if (isstruct(B))
            C=A;
            for (fieldName=fieldnames(B).')
                fieldName=fieldName{1};
                fieldB=B.(fieldName);
                if (isfield(C,fieldName))
                    fieldC=structUnion(C.(fieldName),fieldB);
                    C.(fieldName)=fieldC;
                else
                    C.(fieldName)=fieldB;
                end
            end
        else
            error('incompatible with structure');
        end
    elseif (iscell(A))
        if (iscell(B))
            C=A;
            for element=B
                if (isfield(C,fieldName))
                    fieldC=structUnion(C.(fieldName),fieldB);
                    C.(fieldName)=fieldC;
                else
                    C.(fieldName)=fieldB;
                end
            end
        else
            error('incompatible with cell array');
        end
    elseif (isnumeric(A) || ischar(A))
        if ((isnumeric(A)&&isnumeric(B)) || (isscalar(A)&&isscalar(B)))
            C=B;
        else
            error('incompatible leaf element');
        end
    end
end