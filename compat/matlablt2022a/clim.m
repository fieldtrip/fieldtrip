function clim(varargin)

% This is an drop-in replacement version of the MathWorks function CLIM
% which has been introduced in MATLAB 2022a. This replacement function
% allows CLIM also to be called on older MATLAB versions.
%
% The directory containing this function shoudl only be added to the path
% of MATLAB versions prior to 2022a.

% Copyright (C) 2024, Robert Oostenveld

caxis(varargin{:});
