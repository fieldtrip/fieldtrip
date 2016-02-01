function [params, result] = ortho_default(params, W, w)
%  DSS No orthogonalization function
%  Does mere normalisation
%   W = orthof(W)     for symmetric dss
%   w = orthof(W, w)  for deflation dss
%     W Matrix with projection vectors as rows. For deflation
%       algorithm only previously calculated projections are given.
%     w Currently iterated projection. Only for deflation algorithm.

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Licensed under the Creative Commons Attribution-NonCommercial-ShareAlike
% License. http://creativecommons.org/licenses/by-nc-sa/2.0/.
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<2
    params.name = 'No orthogonalization';
    params.description = 'Description of this function.';
    return;
end

disp('no orthogonalisation!');
if nargin>2
  % per component orthogonalization
  result = w / norm(w);
else
  % symmetric orthogonalization  
  result = diag(1./ sqrt(diag(W * W'))) * W;
end
