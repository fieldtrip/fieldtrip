% NANMIN provides a replacement for MATLAB's nanmin that is almost
% compatible.

function [varargout] = nanmin(varargin)

% the standard MATLAB min function ignores NaNs from at least version 7.0.4 onwards
[varargout{1:nargout}] = min(varargin{:});
