function writematrix(datamatrix, filename, varargin)

% This is a very limited overloaded version of the MathWorks function WRITEMATRIX 
% which has been introduced in MATLAB 2019a. This version does not offer full 
% backwards compatibility, it is only capable of writing text files using 
% the (now deprecated) MathWorks function DLMWRITE under the hood. 
%
% Use as:
%   writematrix(datamatrix, filename, varargin)
% 
% Please read the code to see which 'varargin' arguments are functional

% Copyright (C) 2022 Jan-Mathijs Schoffelen

sel = find(strcmpi(varargin, 'FileType'));
if ~isempty(sel)
  % remove the new-style argument
  varargin([sel sel+1]) = [];
end

sel = find(strcmpi(varargin, 'WriteMode'));
if ~isempty(sel)
  % remove the new-style argument
  varargin([sel sel+1]) = [];
  % prepend the old-style argument
  varargin = ['-append', varargin];
end

dlmwrite(filename, datamatrix, varargin{:})
