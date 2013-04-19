function bbp = bbprctile(x,p,B,w)
% BBPRCTILE  Bayesian bootstrap percentile
%
%    Description:
%    bbp = bbprctile(x,p,B,w)
%    x   = Mx1 matrix of data
%    p   = Px1 or 1xP vector of percentiles
%    B   = number of bootstrap replicates
%    w   = Mx1 vector of weights 

% Copyright (c) Aki Vehtari, 1998-2004
%
% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

error('No mex-file for this architecture. See Matlab help and convert.m in ./linuxCsource or ./winCsource for help.')
