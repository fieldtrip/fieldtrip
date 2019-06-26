function h = uimage(varargin)
%UIMAGE  Display image with uneven axis.
%   UIMAGE(X,Y,C) displays matrix C as an image, using the vectors X and
%   Y to specify the X and Y coordinates. X and Y may be unevenly spaced
%   vectors, but must be increasing. The size of C must be LENGTH(Y)*
%   LENGTH(X). (Most probably you'll want to display C' instead of C).
%
%   Contrary to Matlab's original IMAGE function, here the vectors X and Y
%   do not need to be linearly spaced. Whereas IMAGE linearly interpolates
%   the X-axis between X(1) and X(end), ignoring all other values (idem
%   for Y), UIMAGE allows for X and/or Y to be unevenly spaced vectors, by
%   locally stretching the matrix C (ie, by duplicating some elements of C)
%   for larger X and/or Y intervals.
%
%   The syntax for UIMAGE(X,Y,C,...) is the same as IMAGE(X,Y,C,...)
%   (all the remaining arguments, eg 'PropertyName'-PropertyValue pairs,
%   are passed to IMAGE). See IMAGE for details.
%
%   Use UIMAGESC to scale the data using the full colormap. The syntax for
%   UIMAGESC(X,Y,C,...) is the same as IMAGESC(X,Y,C,...).
%
%   Typical uses:
%      - Plotting a spatio-temporal diagram (T,X), with unevenly spaced
%      time intervals for T (eg, when some values are missing, or when
%      using a non-constant sampling rate).
%      - Plotting a set of power spectra with frequency in log-scale.
%
%   h = UIMAGE(X,Y,C,...) returns a handle to the image.
%
%   Example:
%     c = randn(50,20);         % Random 50x20 matrix
%     x = logspace(1,3,50);     % log-spaced X-axis, between 10 and 1000
%     y = linspace(3,8,20);     % lin-spaced Y-axis, between 3 and 8
%     uimagesc(x,y,c');         % displays the matrix
%
%   F. Moisy
%   Revision: 1.03,  Date: 2006/06/14.
%
%   See also IMAGE, IMAGESC, UIMAGESC.
%
% This function is downloaded on Oct 24th 2008 from www.mathworks.com/matlabcentral/fileexchange/11368

% History:
%   2006/06/12: v1.00, first version.
%   2006/06/14: v1.03, minor bug fixed; works in ML6.
%
% $Id$

narginchk(3,inf);

% maximum number of matrix elements to interpolate the uneven axis
% (typically between 500 and 5000):
nmax = 2000; 

x = varargin{1};
y = varargin{2};
c = varargin{3};

if any(diff(x)<=0) || any(diff(y)<=0)
    error('The X and Y axis should be increasing.');
end

dx = min(diff(x));                   % smallest interval for X
dy = min(diff(y));                   % smallest interval for Y

% test if X and Y are linearly spaced (to within 10^-12):
evenx = all(abs(diff(x)/dx-1)<1e-12);     % true if X is linearly spaced
eveny = all(abs(diff(y)/dy-1)<1e-12);     % true if Y is linearly spaced


if evenx && eveny         % X and Y both evenly spaced

    xe = x;
    ye = y;
    ce = c;

elseif evenx && ~eveny    % X even and Y uneven
    
    nx = length(x);
    xe = x;

    ny = ceil(1 + (y(end) - y(1))/dy);   % number of points for Y
    ny = min(ny, nmax);
    ye = linspace(y(1), y(end), ny);

    ce = zeros(ny,nx);

    for j=1:ny
        indj = find(y<=ye(j));
        ce(j,1:nx) = c(indj(end), 1:nx);
    end;

elseif ~evenx && eveny    % X uneven and Y even
    
    nx = ceil(1 + (x(end) - x(1))/dx);   % number of points for X
    nx = min(nx, nmax);
    xe = linspace(x(1), x(end), nx);

    ny = length(y);
    ye = y;

    ce = zeros(ny,nx);

    for i=1:nx
        indi = find(x<=xe(i));
        ce(1:ny,i) = c(1:ny, indi(end));
    end;

elseif ~evenx && ~eveny   % X and Y both uneven
    
    nx = ceil(1 + (x(end) - x(1))/dx);   % number of points for X
    nx = min(nx, nmax);
    xe = linspace(x(1), x(end), nx);

    ny = ceil(1 + (y(end) - y(1))/dy);   % number of points for Y
    ny = min(ny, nmax);
    ye = linspace(y(1), y(end), ny);

    ce = zeros(ny,nx);
    
    for i=1:nx
        for j=1:ny
            indi = find(x<=xe(i));
            indj = find(y<=ye(j));
            ce(j,i) = c(indj(end),indi(end));
        end;
    end;

end

hh = image(xe, ye, ce, varargin{4:end});

if nargout>0
    h = hh;
end
