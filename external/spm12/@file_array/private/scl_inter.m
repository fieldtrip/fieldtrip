function varargout = scl_inter(varargin)
% file_array's scl_inter property
% For getting the value
% dat = scl_inter(obj)
%
% For setting the value
% obj = scl_inter(obj,dat)
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
dat = obj.scl_inter;


%==========================================================================
% function obj = asgn(obj,dat)
%==========================================================================
function obj = asgn(obj,dat)
if isnumeric(dat) % && numel(dat)<=1,
    obj.scl_inter = double(dat);
else
    error('"scl_inter" must be numeric.');
end