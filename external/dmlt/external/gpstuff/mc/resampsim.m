function s = resampsim(p,m,n)
%RESAMPSIM Simple random resampling
%
%   Description:
%   S = RESAMPSIM(P) returns a new set of indices according to 
%   the probabilities P. P is array of probabilities, which are
%   not necessarily normalized, though they must be non-negative,
%   and not all zero. The size of S is the size of P.
%
%   S = RESAMPSIM(P,M,N) returns M by N matrix.
%
%   S = RESAMPSIM(P,M) returns M by M matrix.
%
%   Note that residual, stratified and deterministic resampling all
%   have smaller variance.
%
%   Simple random resampling samples indices randomly according
%   to the probabilities P. See, e.g., Liu, J. S., Monte Carlo
%   Strategies in Scientific Computing, Springer, 2001, p. 72.
%
%   See also RESAMPRES, RESAMPSTR, RESAMPDET

% Copyright (c) 2003-2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin<2
    [m,n] = size(p);
elseif nargin==2
    n = m;
end
pc=cumsum(p(:));
pc=pc./pc(end);
s=binsgeq(pc,rand([m,n]));
