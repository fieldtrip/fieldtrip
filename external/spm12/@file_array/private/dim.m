function varargout = dim(varargin)
% file_array's dimension property
% For getting the value
% dat = dim(obj)
%
% For setting the value
% obj = dim(obj,dat)
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
dat = obj.dim;


%==========================================================================
% function obj = asgn(obj,dat)
%==========================================================================
function obj = asgn(obj,dat)
if isnumeric(dat) && all(dat>=0) && all(rem(dat,1)==0)
    dat = [double(dat(:)') 1 1];
    lim = max([2 find(dat~=1)]);
    dat = dat(1:lim);
    obj.dim = dat;
else
    error('"dim" must be a vector of positive integers.');
end
