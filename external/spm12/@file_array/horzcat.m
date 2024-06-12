function o = horzcat(varargin)
% Horizontal concatenation of file_array objects
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


o = cat(2,varargin{:});
