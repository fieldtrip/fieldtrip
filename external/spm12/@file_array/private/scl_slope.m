function varargout = scl_slope(varargin)
% file_array's scl_slope property
% For getting the value
% dat = scl_slope(obj)
%
% For setting the value
% obj = scl_slope(obj,dat)
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
dat = obj.scl_slope;


%==========================================================================
% function obj = asgn(obj,dat)
%==========================================================================
function obj = asgn(obj,dat)
if isnumeric(dat) % && numel(dat)<=1,
    obj.scl_slope = double(dat);
else
    error('"scl_slope" must be numeric.');
end
