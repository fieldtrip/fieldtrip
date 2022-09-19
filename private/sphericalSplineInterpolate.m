function [W,Gss,Gds,Hds]=sphericalSplineInterpolate(src,dest,lambda,order,type,tol)
%interpolate matrix for spherical interpolation
%
% W = sphericalSplineInterpolate(src,dest,lambda,order,type,tol)
%
% Inputs:
%  src    - [3 x N] old electrode positions
%  dest   - [3 x M] new electrode positions
%  lambda - [float] regularisation parameter for smoothing the estimates (1e-5)
%  order  - [float] order of the polynomial interplotation to use (4)
%  type - [str] one of;                                         ('spline')
%             'spline' - spherical Spline 
%             'slap'   - surface Laplician (aka. CSD)
%  tol    - [float] tolerance for the legendre poly approx        (1e-7)
% Outputs:
%  W      - [M x N] linear mapping matrix between old and new co-ords
%
% Based upon the paper: Perrin89

% Copyright 2009-     by Jason D.R. Farquhar (jdrf@zepler.org)
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents, provided this copyright
% notice is retained, and note is made of any changes that have been
% made. This software and documents are distributed without any
% warranty, express or implied. 

if ( nargin < 3 || isempty(lambda)) lambda=1e-5; end
if ( nargin < 4 || isempty(order)) order=4; end
if ( nargin < 5 || isempty(type)) type='spline'; end
if ( nargin < 6 || isempty(tol)) tol=eps; end

% map the positions onto the sphere (not using repop, by JMH)
src  = src ./repmat(sqrt(sum(src.^2)),  size(src, 1), 1);  % src   = repop(src,'./',sqrt(sum(src.^2)));
dest = dest./repmat(sqrt(sum(dest.^2)), size(dest, 1), 1); % dest  = repop(dest,'./',sqrt(sum(dest.^2)));

%calculate the cosine of the angle between the new and old electrodes. If
%the vectors are on top of each other, the result is 1, if they are
%pointing the other way, the result is -1
cosSS = src'*src;  % angles between source positions
cosDS = dest'*src; % angles between destination positions

% Compute the interpolation matrix to tolerance tol
[Gss]      = interpMx(cosSS,order,tol);  % [nSrc x nSrc]
[Gds, Hds] = interpMx(cosDS,order,tol);  % [nDest x nSrc]

% Include the regularisation
if ( lambda>0 ) Gss = Gss+lambda*eye(size(Gss)); end

% Compute the mapping to the polynomial coefficients space % [nSrc+1 x nSrc+1]
% N.B. this can be numerically unstable so use the PINV to solve..
muGss=1; %median(diag(Gss)); % used to improve condition number when inverting. Probably uncessary
%C = [      Gss            muGss*ones(size(Gss,1),1)];
C = [      Gss            muGss*ones(size(Gss,1),1);...
      muGss*ones(1,size(Gss,2))       0];
iC = pinv(C);

% Compute the mapping from source measurements and positions to destination positions
if ( strcmpi(type,'spline') )
  W = [Gds ones(size(Gds,1),1).*muGss]*iC(:,1:end-1); % [nDest x nSrc]
elseif (strcmpi(type,'slap'))
  W = Hds*iC(1:end-1,1:end-1); %(:,1:end-1); % [nDest x nSrc]
end
return;
%--------------------------------------------------------------------------
function [G,H]=interpMx(cosEE,order,tol)
% compute the interpolation matrix for this set of point pairs
if ( nargin < 3 || isempty(tol) ) tol=1e-10; end
G=zeros(size(cosEE)); H=zeros(size(cosEE));
for i=1:numel(cosEE)
   x = cosEE(i);
   n=1; Pns1=1; Pn=x;           % seeds for the legendre ploy recurence
   tmp  = ( (2*n+1) * Pn ) / ((n*n+n).^order);
   G(i) = tmp ;         % 1st element in the sum
   H(i) = (n*n+n)*tmp;  % 1st element in the sum
   oGi=inf; dG=abs(G(i)); oHi=inf; dH=abs(H(i));
   for n=2:500 % do the sum
      Pns2=Pns1; Pns1=Pn; Pn=((2*n-1)*x*Pns1 - (n-1)*Pns2)./n; % legendre poly recurance
      oGi=G(i);  oHi=H(i);
      tmp  = ((2*n+1) * Pn) / ((n*n+n).^order);
      G(i) = G(i) + tmp;                     % update function estimate, spline interp     
      H(i) = H(i) + (n*n+n)*tmp;             % update function estimate, SLAP
      dG   = (abs(oGi-G(i))+dG)/2; dH=(abs(oHi-H(i))+dH)/2; % moving ave gradient est for convergence
      %fprintf('%d) dG =%g \t dH = %g\n',n,dG,dH); %abs(oGi-G(i)),abs(oHi-H(i)));
      if ( dG<tol && dH<tol ) break; end           % stop when tol reached
   end
end
G= G./(4*pi);
H= H./(4*pi);
return;
%--------------------------------------------------------------------------
function testCase()
src=randn(3,100); src(3,:)=abs(src(3,:)); src=repop(src,'./',sqrt(sum(src.^2))); % source points on sphere
dest=rand(3,30); dest(3,:)=abs(dest(3,:)); dest=repop(dest,'./',sqrt(sum(dest.^2))); % dest points on sphere
clf;scatPlot(src,'b.');

W=sphericalSplineInterpolate(src,dest);

W=sphericalSplineInterpolate(src(:,1:70),src);
clf;imagesc(W);
clf;jplot(src(1:2,1:70),W(70:75,:)');

z=jf_load('eeg/vgrid/nips2007/1-rect230ms','jh','flip_rc_sep');
z=jf_retain(z,'dim','ch','idx',[z.di(1).extra.iseeg]);
lambda=0; order=4;
chPos=[z.di(1).extra.pos3d];
incIdx=1:size(dest,2)-3; exIdx=setdiff(1:size(dest,2),incIdx); % included + excluded channels
W=sphericalSplineInterpolate(chPos(:,incIdx),chPos,lambda,order); % estimate the removed channels
clf;imagesc(W);
clf;jplot(chPos(:,incIdx),W(exIdx,:)');
% compare estimated with true
X=z.X(:,:,2);
Xest=W*z.X(incIdx,:,2);
clf;mcplot([X(exIdx(1),:);Xest(exIdx(1),:)]')
clf;subplot(211);mcplot(X');subplot(212);mcplot(Xest');

