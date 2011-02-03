function varargout = fname(varargin)
% Format
% For getting the value
% dat = fname(obj)
%
% For setting the value
% obj = fname(obj,dat)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id$



if nargin==2,
    varargout{1} = asgn(varargin{:});
elseif nargin==1,
    varargout{1} = ref(varargin{:});
else
    error('Wrong number of arguments.');
end;
return;

function dat = ref(obj)
dat = obj.fname;
return;

function obj = asgn(obj,dat)
if ischar(dat)
    obj.fname = deblank(dat(:)');
else
    error('"fname" must be a character string.');
end;
return;
