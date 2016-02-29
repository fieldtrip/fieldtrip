function source = smooth_source(source, varargin)

%[SOURCE] = SMOOTH(SOURCE, VARARGIN)
%
% computes location specific 3D gaussian kernels based on a FWHM estimate
%  source should contain the fields 
%    fwhm, specifying for each voxel the FWHM of the smoothing kernel in the xyz-direction
%    pos,  allowing for the units to be correct
%
%  key-value pairs should contain
%    parameter = string, field to be used for the smoothing
%    maxdist   = scalar, maximum distance for filter kernel

source    = ft_datatype_source(source);

parameter = ft_getopt(varargin, 'parameter', 'all');
maxdist   = ft_getopt(varargin, 'maxdist', []);

if ischar(parameter) && strcmp(parameter, 'all')
  error('not yet implemented');	
elseif ischar(parameter)
	parameter = {parameter};
end

if isempty(maxdist)
	maxdist = inf;
end

if islogical(source.inside)
	inside = find(source.inside);
else
	inside = source.inside;
end
ninside = numel(inside);

if isfield(source, 'fwhm')
  fwhm = source.fwhm(inside,:,:);
else
	error('the input data should contain an ''fwhm'' field');
end

dat = cell(1,numel(parameter));
for k = 1:numel(parameter)
  dat{k} = getsubfield(source, parameter);
  dat{k} = dat{k}(inside,:);
end

pos = source.pos(inside,:);

sigma = zeros(ninside, 3, 3);
if numel(fwhm)==ninside
	sigma(:,1,1) = (fwhm./(2.*sqrt(2.*log(2)))).^2;
  sigma(:,2,2) = (fwhm./(2.*sqrt(2.*log(2)))).^2;
  sigma(:,3,3) = (fwhm./(2.*sqrt(2.*log(2)))).^2;
else
	% this should account for more fancy kernels
	error('not yet implemented');
end

tmp = cell(1,numel(parameter));
for m = 1:numel(tmp)
	tmp{m} = zeros(size(dat{k}));
end
for k = 1:ninside
  thispos = pos(k,:);
	dpos    = pos - thispos(ones(ninside,1),:);
	sel     = sqrt(sum(dpos.^2,2))<=maxdist;
	krn     = gaussian3d(squeeze(sigma(k,:,:)), dpos(sel,:));
  krn     = krn./sum(krn);
	for m = 1:numel(dat)
		tmp{m}(k,:) = krn'*dat{m}(sel,:);
	end
end

smoothdat = zeros(size(source.pos,1),size(tmp{1},2));
for m = 1:numel(tmp)
  smoothdat(inside,:) = tmp{m};
  source = setsubfield(source, [parameter{m},'smooth'], smoothdat);
end

%---------------------------------------
function [output] = gaussian3d(C, x, y, z)

% [OUTPUT] = GAUSSIAN3D(C, x, y, z) computes a 3D gaussian volume, as
% specified by the covariance matrix C. x,y,z specify the x, y, and z
% distance to the center of the gaussian

if nargin==2 && size(x,2)==3,
	do3Dvol = false;
elseif nargin==2
	y = x;
	z = x;
	do3Dvol = true;
else
	do3Dvol = true;
end

if do3Dvol
	lx = numel(x);
	ly = numel(y);
	lz = numel(z);
	
	x = reshape(x, [lx 1 1]);
	y = reshape(y, [1 ly 1]);
	z = reshape(z, [1 1 lz]);
	
	rx = ones(1,lx);
	ry = ones(1,ly);
	rz = ones(1,lz);
	
	xy = x(:,ry,:).*y(rx,:,:);
	xz = x(:,:,rz).*z(rx,:,:);
	yz = y(:,:,rz).*z(:,ry,:);
	
	xsq = x.^2;
	ysq = y.^2;
	zsq = z.^2;
	
	%regularize a bit if necessary
	ev = eig(C);
	if ev(3)/ev(1)>1e12,
		s = svd(C);
		C = C+eye(size(C)).*0.001.*s(1);
	else
	end
	
	detC = det(C);
	invC = inv(C);
	A    = 1/((2*pi)^1.5*sqrt(detC));
	
	c11 = -0.5.*invC(1,1).*xsq;
	c12 = -invC(1,2).*xy;
	c13 = -invC(1,3).*xz;
	c22 = -0.5.*invC(2,2).*ysq;
	c23 = -invC(2,3).*yz;
	c33 = -0.5.*invC(3,3).*zsq;
	
	%output = A*exp(-0.5*(repmat(c11.*x.^2, [1 ly lz]) + ...
	%                     repmat(c22.*y.^2, [lx 1 ly]) + ...
	%		     repmat(c33.*z.^2, [lx ly 1]) + ...
	%		     repmat(2.*c12.*xy, [1 1 lz]) + ...
	%		     repmat(2.*c13.*xz, [1 ly 1]) + ...
	%		     repmat(2.*c23.*yz, [lx 1 1])));
	%output = A*exp(-0.5*(c11.*xsq(:,ones(1,ly),ones(1,lz)) + ...
	%                     c22.*ysq(ones(1,lx),:,ones(1,lz)) + ...
	%		     c33.*zsq(ones(1,lx),ones(1,ly),:) + ...
	%		     2.*c12.*xy(:,:,ones(1,lz)) + ...
	%		     2.*c13.*xz(:,ones(1,ly),:) + ...
	%		     2.*c23.*yz(ones(1,lx),:,:)));
	
	%output = A*exp((c11(:,ry,rz) + ...
	%                c22(rx,:,rz) + ...
	%	        c33(rx,ry,:) + ...
	%	        c12(:,:,rz)  + ...
	%	        c13(:,ry,:)  + ...
	%	        c23(rx,:,:)));
	
	%c12 = c12+c11(:,ry);
	%c23 = c23+c22(:,:,rz);
	%c13 = c13+c33(rx,:,:);
	%
	%output = A*exp( c12(:,:,rz)  + ...
	%	        c13(:,ry,:)  + ...
	%	        c23(rx,:,:));
	
	c12 = exp(c12+c11(:,ry));
	c23 = exp(c23+c22(:,:,rz));
	c13 = exp(c13+c33(rx,:,:));
	
	output = A*c12(:,:,rz).*c13(:,ry,:).*c23(rx,:,:);
else
	% treat as a Nx3 list of x,y,z distances, and output Nx1 kernel
	%regularize a bit if necessary
	ev = eig(C);
	if ev(3)/ev(1)>1e12,
		s = svd(C);
		C = C+eye(size(C)).*0.001.*s(1);
	else
	end
	
	detC = det(C);
	invC = inv(C);
	A    = 1/((2*pi)^1.5*sqrt(detC));
	
	% efficiently compute x*invC*x'
	sandwich = x(:,1).^2.*invC(1,1)+x(:,2).^2.*invC(2,2)+x(:,3).^2.*invC(3,3)+...
		         x(:,1).*x(:,2).*invC(1,2).*2 + ...
						 x(:,1).*x(:,3).*invC(1,3).*2 + ...
						 x(:,2).*x(:,3).*invC(2,3).*2;
					 
	
	output = A.*exp(-0.5.*sandwich);
end
