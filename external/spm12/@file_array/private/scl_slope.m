function varargout = scl_slope(varargin)
% Format
% For getting the value
% dat = scl_slope(obj)
%
% For setting the value
% obj = scl_slope(obj,dat)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: scl_slope.m 1143 2008-02-07 19:33:33Z spm $



if nargin==2,
    varargout{1} = asgn(varargin{:});
elseif nargin==1,
    varargout{1} = ref(varargin{:});
else
    error('Wrong number of arguments.');
end;
return;

function dat = ref(obj)
dat = obj.scl_slope;
return;

function obj = asgn(obj,dat)
if isnumeric(dat), % && numel(dat)<=1,
    obj.scl_slope = double(dat);
else
    error('"scl_slope" must be numeric.');
end;
return;
