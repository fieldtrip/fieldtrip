function o = vertcat(varargin)
% Vertical concatenation of file_array objects.
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


o = cat(1,varargin{:});
