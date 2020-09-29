function objnew = CopyHandles(varargin)
% Syntax:
%       objnew = CopyHandles(obj)
%       objnew = CopyHandles(obj, objnew)
%
% Generate copy of array of handle objects from obj to objnew
%
%

if nargin==1
    obj = varargin{1};
elseif nargin==2
    obj    = varargin{1};
    objnew = varargin{2};
end

if isempty(obj)
    objnew = obj;
    return
end

if isa(obj, 'handle')
    % Generate class instance with shallow copy
    objnew = obj(1).copy;
    for ii=2:length(obj)
        objnew(ii) = obj(ii).copy;
    end
else
    objnew = obj;
end
    
props = properties(obj(1));

% Deep copy properties from obj to objnew
for ii=1:length(obj)   
    for iProp=1:length(props)
        if isa( eval( sprintf('obj(ii).%s', props{iProp}) ), 'handle')
            for kk=1:eval( sprintf('length(obj(ii).%s)', props{iProp}) )
                eval( sprintf('objnew(ii).%s(kk) = CopyHandles(obj(ii).%s(kk));', props{iProp}, props{iProp}) );
            end
        else
            % Non-handle object assignment does copy-by-value which is deep copy by definition 
            eval( sprintf('objnew(ii).%s = obj(ii).%s;', props{iProp}, props{iProp}) );
        end
    end
end

