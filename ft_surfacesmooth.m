function bnd = ft_surfacesmooth(cfg,bnd)
% FT_SURFACESMOOTH returns a smoother version of the surface input
% 
% Use as
%   cfg = [];
%   cfg.method = 'spharm';
%   bnd = ft_surfacesmooth(cfg,bnd)
% 
% See also FT_SURFACECHECK, FT_SURFACEEXTRACT

% Copyright (C) 2011, Cristiano Micheli, Robert Oostenveld
% 
% $Id: $

ft_defaults

method  = ft_getopt(cfg, 'method', 'spharm');
sphord  = ft_getopt(cfg, 'sphord', 12); % sphord=12, see van 't Ent 

switch method
  
  case 'spharm'
    cfg = [];
    cfg.isclosed = 'yes';
    ft_surfacecheck(cfg,bnd);
    [bnd] = spherical_harmonic_mesh(bnd, sphord);

  otherwise
    error('unsupported method "%s"', cfg.method); 
    
end

% Attention .tri is flipped here (see end), otherwise normals' directions change
function [bnd] = spherical_harmonic_mesh(bnd, nl)
% SPHERICAL_HARMONIC_MESH realizes the smoothed version of a mesh contained in the first
% argument bnd. The boundary argument (bnd) contains typically 2 fields
% called .pnt and .tri referring to vertices and triangulation of a mesh.
% The degree of smoothing is given by the order in the second input: nl
%
% Use as
%   bnd2 = spherical_harmonic_mesh(bnd, nl);

bnd2 = [];

% calculate midpoint
Or = mean(bnd.pnt);
% rescale all points
x = bnd.pnt(:,1) - Or(1);
y = bnd.pnt(:,2) - Or(2);
z = bnd.pnt(:,3) - Or(3);
X = [x(:), y(:), z(:)];

% convert points to parameters
% [T,P] = points2param(X);
x = X(:,1); y = X(:,2); z = X(:,3);
[phi,theta,R] = cart2sph(x,y,z);
theta = theta + pi/2;
phi   = phi + pi;
T     = theta(:);
P     = phi(:);

% basis function
B     = shlib_B(nl, T, P);
Y     = B'*X;

% solve the linear system
a     = pinv(B'*B)*Y;

% build the surface
Xs = zeros(size(X));
for l = 0:nl-1
  for m = -l:l
    if m<0
      Yml = shlib_Yml(l, abs(m), T, P);
      Yml = (-1)^m*conj(Yml);
    else
      Yml = shlib_Yml(l, m, T, P);
    end
    indx = l^2 + l + m+1;
    Xs = Xs + Yml*a(indx,:);
  end
  fprintf('%d of %d\n', l, nl);
end
% Just take the real part
Xs = real(Xs);
% reconstruct the matrices for plotting.
xs = reshape(Xs(:,1), size(x,1), size(x,2));
ys = reshape(Xs(:,2), size(x,1), size(x,2));
zs = reshape(Xs(:,3), size(x,1), size(x,2));

bnd2.pnt = [xs(:)+Or(1) ys(:)+Or(2) zs(:)+Or(3)];
[bnd2.tri] = projecttri(bnd.pnt);
bnd2.tri = fliplr(bnd2.tri);

function B = shlib_B(nl, theta, phi)
% function B = shlib_B(nl, theta, phi)
%
% Constructs the matrix of basis functions
% where b(i,l^2 + l + m) = Yml(theta, phi)
%
% See also: shlib_Yml.m, shlib_decomp.m, shlib_gen_shfnc.m
%
% Dr. A. I. Hanna (2006).
B = zeros(length(theta), nl^2+2*nl+1);
for l = 0:nl-1
  for m = -l:l
    if m<0
      Yml = shlib_Yml(l, abs(m), theta, phi);
      Yml = (-1)^m*conj(Yml);
    else
      Yml = shlib_Yml(l, m, theta, phi);
    end
    indx = l^2 + l + m+1;
    B(:, indx) = Yml;
  end
end
return;

function Yml = shlib_Yml(l, m, theta, phi)
% function Yml = shlib_Yml(l, m, theta, phi)
%
% A matlab function that takes a given order and degree, and the matrix of
% theta and phi and constructs a spherical harmonic from these. The
% analogue in the 1D case would be to give a particular frequency.
%
% Inputs:
%  m - order of the spherical harmonic
%  l - degree of the spherical harmonic
%  theta - matrix of polar coordinates \theta \in [0, \pi]
%  phi - matrix of azimuthal cooridinates \phi \in [0, 2\pi)
%
% Example:
%
% [x, y] = meshgrid(-1:.1:1);
% z = x.^2 + y.^2;
% [phi,theta,R] = cart2sph(x,y,z);
% Yml = spharm_aih(2,2, theta(:), phi(:));
%
% See also: shlib_B.m, shlib_decomp.m, shlib_gen_shfnc.m
%
% Dr. A. I. Hanna (2006)
Pml=legendre(l,cos(theta));
if l~=0
  Pml=squeeze(Pml(m+1,:,:));
end
Pml = Pml(:);
% Yml = sqrt(((2*l+1)/4).*(factorial(l-m)/factorial(l+m))).*Pml.*exp(sqrt(-1).*m.*phi);
% new:
Yml = sqrt(((2*l+1)/(4*pi)).*(factorial(l-m)/factorial(l+m))).*Pml.*exp(sqrt(-1).*m.*phi);
return;
