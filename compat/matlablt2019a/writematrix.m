function writematrix(datamatrix, filename, varargin)

% This is a extremely limited backward compatible overloaded version of
% writematrix with extremely limited functionality. The writematrix function
% has been introduced in MATLAB 2019a. This version is only capable of 
% % writing text files, using dlmwrite
% under the hood. dlmwrite is not recommended by Mathworks, where the
% intended replacement function is writematrix
%
% Use as:
%   writematrix(datamatrix, filename, varargin)
% 
% Please read the code to see which 'varargin' arguments are functional

% Copyright (C) 2022 Jan-Mathijs Schoffelen

sel = find(strcmpi(varargin, 'filetype'));
if isempty(sel)
  filetype = 'text';
else
  filetype = varargin{sel+1};
  varargin([sel sel+1]) = [];
end
dlmwrite(filename, datamatrix, varargin{:})
