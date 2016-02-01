function R = gamrand1(A, B)
%GAMRAND1 Random matrices from gamma distribution.
%
%   R = GAMRAND1(A,B) returns a matrix of random numbers chosen   
%   from the gamma distribution with parameters A and B.
%   Both parameters have to be scalar. 
% 
%   Note: Parameterization as in (Neal, 1996).
%      A is mean of the distribution
%      B is degrees of freedom
%
%	See also INVGAMRAND

% Copyright (c) 1999 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

error('No mex-file for this architecture. See Matlab help and convert.m in ./linuxCsource or ./winCsource for help.')
