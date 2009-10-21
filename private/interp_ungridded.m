function varargout = interp_ungridded(varargin);

% INTERP_UNGRIDDED computes an interpolation matrix for two ungridded
% clouds of 3-D points
%
% Use as
%   [val] = interp_ungridded(pnt1, val1, pnt2, ...)
% where
%   pnt1    Nx3 matrix with the vertex positions
%   val1    Nx1 vector with the values on each vetrex (can also be multiple columns)
%   pnt2    Mx3 matrix with the vertex positions onto which the data should
%           be interpolated
%
% Alternatively to get the interpolation matrix itself, you can use it as
%   [interpmat, distmat] = interp_ungridded(pnt1, pnt2, ...)
% 
%
% Optional arguments are specified in key-value pairs and can be
%    projmethod   = 'nearest', 'sphere_avg', 'sphere_weighteddistance'
%    sphereradius = number
%    distmat      = NxM matrix with precomputed distances
%    tri1         = triangulation for the first set of vertices

% Copyright (C) 2007, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% $Log: interp_ungridded.m,v $
% Revision 1.1  2007/02/07 07:42:06  roboos
% new implementation to be used in teh new sourceplot, mainly based on code from the old surfaceplot
%

if nargin<3
  error('Not enough input arguments.');
end

% get the fixed input arguments
if isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3})
  pnt1 = varargin{1};
  val1 = varargin{2};
  pnt2 = varargin{3};
  varargin = varargin(4:end);
  hasval = 1;
elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) && ischar(varargin{3})
  pnt1 = varargin{1};
  pnt2 = varargin{2};
  varargin = varargin(3:end);
  hasval = 0;
end

% get the optional arguments
projmethod   = keyval('projmethod',    varargin);   % required
sphereradius = keyval('sphereradius',  varargin);   % required for some projection methods
distmat      = keyval('distmat',       varargin);   % will be computed if not present
tri1         = keyval('tri1', varargin);            % not yet implemented

npnt1 = size(pnt1, 1);
npnt2 = size(pnt2, 1);

if isempty(distmat)
  %------compute a distance matrix
  switch projmethod
    case 'nearest'
      if ~isempty(sphereradius)
        warning('sphereradius is not used for projmethod''nearest''');
      end
      % determine the nearest voxel for each surface point
      ind = find_nearest(pnt2, pnt1, 5);
      distmat = sparse(1:npnt2, ind, ones(size(ind)), npnt2, npnt1);

    case {'sphere_avg', 'sphere_weighteddistance', 'sphere_weightedprojection'}
      if isempty(sphereradius)
        error('sphereradius should be specified');
      end
      % compute the distance between voxels and each surface point
      dpnt1sq  = sum(pnt1.^2,2); % squared distance to origin
      dpnt2sq  = sum(pnt2.^2,2); % squared distance to origin
      maxnpnt = double(npnt2*ceil(4/3*pi*(sphereradius/max(dimres))^3)); % initial estimate of nonzero entries
      distmat = spalloc(npnt2, npnt1, maxnpnt);
      progress('init', 'textbar', 'computing distance matrix');
      for j = 1:npnt2
        progress(j/npnt2);
        d   = sqrt(dpnt2sq(j) + dpnt1sq - 2 * pnt1 * pnt2(j,:)');
        sel = find(d<sphereradius);
        distmat(j, sel) = single(d(sel)) + eps('single');
      end
      progress('close');

    otherwise
      error('unsupported projection method');
  end % case projmethod
end % if isempty distmat

%------do something with the distance matrix
switch projmethod
  case 'nearest'
    projmat         = distmat;

  case 'sphere_avg'
    projmat         = distmat;
    [ind1, ind2, d] = find(projmat);
    nnz             = full(sum(spones(projmat),2));
    for k = 1:length(ind1)
      projmat(ind1(k),ind2(k)) = 1./nnz(ind1(k));
    end

  case 'sphere_weighteddistance'
    projmat         = distmat;
    [ind1, ind2, d] = find(projmat);
    projmat         = sparse(ind1, ind2, 1./d, npnt, npnt1);
    [ind1, ind2, d] = find(projmat);
    normnz          = sqrt(full(sum(projmat.^2, 2)));
    projmat         = sparse(ind1, ind2, d./normnz(ind1), npnt, npnt1);

  case 'sphere_weightedprojection'
    % JM had something in mind for this, but it is not yet implemented
    error('unsupported projection method');

  otherwise
    error('unsupported projection method');
end  % case projmethod

if hasval
  % return the interpolated values
  varargout{1} = projmat * val1;
else
  % return the interpolation and the distance matrix
  varargout{1} = projmat;
  varargout{2} = distmat;
end

