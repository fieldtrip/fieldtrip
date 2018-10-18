% NANMAX provides a replacement for MATLAB's nanmax that is almost
% compatible.

function [varargout] = nanmax(varargin)

% the standard MATLAB max function ignores NaNs from at least version 7.0.4 onwards
[varargout{1:nargout}] = max(varargin{:});
