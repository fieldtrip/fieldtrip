function [cfg] = keyval2cfg(varargin);

% KEYVAL2CFG converts between a structure and a cell-array with key-value
% pairs which can be used for optional input arguments. 
% 
% Use as
%   [cfg] = keyval2cfg(varargin)

if iscell(varargin) && length(varargin)==1
  varargin = varargin{1};
end

% assign the optional key-value arguments to a configuration structure
var = varargin(1:2:length(varargin));   % get the odd arguments
val = varargin(2:2:length(varargin));   % get the even arguments
cfg = cell2struct(val(:), var(:), 1);
