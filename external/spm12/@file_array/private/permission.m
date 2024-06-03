function varargout = permission(varargin)
% file_array's permission property
% For getting the value
% dat = permission(obj)
%
% For setting the value
% obj = permission(obj,dat)
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


if nargin==2
    varargout{1} = asgn(varargin{:});
elseif nargin==1
    varargout{1} = ref(varargin{:});
else
    error('Wrong number of arguments.');
end


%==========================================================================
% function dat = ref(obj)
%==========================================================================
function dat = ref(obj)
dat = obj.permission;


%==========================================================================
% function obj = asgn(obj,dat)
%==========================================================================
function obj = asgn(obj,dat)
if ischar(dat)
    tmp = lower(deblank(dat(:)'));
    switch tmp
        case 'ro'
        case 'rw'
        otherwise
            error('Permission must be either "ro" or "rw".');
    end
    obj.permission = tmp;
else
    error('"permission" must be a character string.');
end
