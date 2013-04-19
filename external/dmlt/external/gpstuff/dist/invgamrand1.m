function r = invgamrand1(a, b, varargin)
%INVGAMRAND1 Random matrices from inverse gamma distribution
%
%   R = INVGAMRAND1(A,B) returns a matrix of random numbers chosen   
%   from the inverse gamma distribution with parameters A and B.
%   The size of R is the common size of A and B if both are matrices.
%   Both parameters have to be a scalar.
% 
%   Note: Parameterization as in (Neal, 1996).
%      A is mean of the distribution
%      B is degrees of freedom
%   
%	See also INVGAMRAND, GAMRAND1

% Copyright (c) 1999 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

error('No mex-file for this archtitecture. See Matlab help and convert.m in ./linuxCsource or ./winCsource for help.')
