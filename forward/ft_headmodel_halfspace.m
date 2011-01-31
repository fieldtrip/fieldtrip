function vol = ft_headmodel_halfspace(geom, Pvoid, varargin)

% FT_HEADMODEL_HALFSPACE creates an EEG volume conduction model that
% is described with an infinite conductive halfspace. You can think
% of this as a plane with on one side a infinite mass of conductive
% material (e.g. water) and on the other side non-conductive material
% (e.g. air).
%   geom.pnt = Nx3 vector specifying N points through which a plane is fitted 
%   Pvoid    = 1x3 vector specifying the spatial position of a point lying in the empty halfspace 
%              (this determines the plane normal's direction)
% Use as
%   vol = ft_headmodel_halfspace(geom, Pvoid, varargin)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   the non-connductive halfspace)
%   conductivity  = number,  conductivity value of the conductive halfspace
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

cond  = keyval('conductivity', varargin); 
if isempty(cond), cond = 1; warning('Unknown conductivity value (set to 1)'), end

% the description of this volume conduction model consists of the
% description of the plane, and a point in the void halfspace
vol = [];

if isfield(geom,'pnt')
  pnt = geom.pnt;
  [err,N,P] = fit_plane(pnt(:,1), pnt(:,2), pnt(:,3));
else
  error('geometry format is incorrect')
end

vol.cond  = cond;
vol.Pvoid = Pvoid;
vol.pnt   = P(:)'; % a point that lies in the plane that separates the conductive tissue from the air
vol.ori   = N(:)'; % a unit vector pointing towards the air
vol.ori   = vol.ori/norm(vol.ori);
vol.type  = 'halfspace';

function [Err,N,P] = fit_plane(XData, YData, ZData)
% Fits a plane through a number of points in 3D cartesian coordinates
% 
% This and next function are wrapper functions to some pieces of the code from 
% the Statistics Toolbox demo titled "Fitting an Orthogonal 
% Regression Using Principal Components Analysis" 
% (http://www.mathworks.com/products/statistics/
%  demos.html?file=/products/demos/shipping/stats/orthoregdemo.html),
% which is Copyright by the MathWorks, Inc.
X(:,1) = XData(:,1);
X(:,2) = YData(:,1);
X(:,3) = ZData(:,1);
meanX  = mean(X,1);
[coeff,score] = princomp(X);
normal = coeff(:,3);
[n,p]  = size(X);
meanX  = mean(X,1);
error  = abs((X - repmat(meanX,n,1))*normal);
% outputs
Err = sum(error);
N   = normal;
P   = meanX;

function [coeff,score]=princomp(x)
% Principal component analysis, rewritten for compatibility
if ~isempty(x)
  [n,p] = size(x);
  % Center X by subtracting off column means
  x0 = bsxfun(@minus,x,mean(x,1));
  r = min(n-1,p); % max possible rank of X0
  [U,sigma,coeff] = svd(x0,0); 
  if n == 1 % sigma might have only 1 row
      sigma = sigma(1);
  else
      sigma = diag(sigma);
  end
  score = bsxfun(@times,U,sigma'); % == x0*coeff
  sigma = sigma ./ sqrt(n-1);
  % When X has at least as many variables as observations, eigenvalues
  % n:p of S are exactly zero.
  if n <= p
    sigma(n:p,1) = 0; % make sure this extends as a column
    score(:,n:p) = 0;
  end
end
