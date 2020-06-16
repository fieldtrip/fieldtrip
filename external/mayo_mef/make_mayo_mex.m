function  make_mayo_mex()
% MAKE_MAYO_MEX make mex binary necessary for reading MEF dataset
% 
% Syntax:
%   make_mayo_mex
% 
% Input(s):
%
% Output(s):
%
% Example:
%
% Note:
%
% References:
%
% See also .

% Copyright 2020 Richard J. Cui. Created: Fri 05/15/2020 10:33:00.474 AM
% $ Revision: 0.1 $  $ Date: Fri 05/15/2020 10:33:00.474 AM $
%
% Multimodel Neuroimaging Lab (Dr. Dora Hermes)
% Mayo Clinic St. Mary Campus
% Rochester, MN 55905, USA
%
% Email: richard.cui@utoronto.ca (permanent), Cui.Jie@mayo.edu (official)

cur_dir = pwd;

mayo_mef = fileparts(mfilename('fullpath')); % directory of make_mayo_mex.m assumed in mayo_mef
cd([mayo_mef,filesep,'mex_mef'])
make_mex_mef

cd(cur_dir)

end

% [EOF]