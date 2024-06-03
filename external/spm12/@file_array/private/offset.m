function varargout = offset(varargin)
% file_array's offset property
% For getting the value
% dat = offset(obj)
%
% For setting the value
% obj = offset(obj,dat)
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
dat = obj.offset;


%==========================================================================
% function obj = asgn(obj,dat)
%==========================================================================
function obj = asgn(obj,dat)
if isnumeric(dat) && numel(dat)==1 && dat>=0 && rem(dat,1)==0
    obj.offset = double(dat);
else
    error('"offset" must be a positive integer.');
end
