function [center,radius]=sphfit(vc,Ni,stopth)
% fits a sphere to a set of surface points
% usage: [center,radius]=sphfit(vc,Ni,stopth)
%
% input:
% vc        nx3 matrix, where each row represents the location
%           of one surface point. vc can have more than 3 columns
%           (e.g. orientations) - then only the first 3 columns are used
% Ni        number of iteration requested, 20 by default
% stopth    stopping threshold used for the relative difference between 2
%           succeessive estimates of the radius. Fixed 10^-6 by default
%           If stopth<0, then no stopping till the end of the Ni iterations
%
% center  1x3 vector denoting the center
% radius  scalar denoting the radius
%
% written by Guido Nolte
% updated by Christophe Phillips, 2009/1/19
%   - add the number of iterations as input, use 5 as default
%   - add a stopping criteria based on the relative difference between 2
%       successive estimates of the radius.
%       If rel_diff<1e-6 (default of use fixed), then break out.

if nargin<2 || isempty(Ni)
  Ni = 20;
end

if nargin<3 || isempty(stopth)
  stopth = 1e-6;
end

vc=vc(:,1:3);
[nvc,ndum]=size(vc);

center_0=mean(vc);
vcx=vc-repmat(center_0,nvc,1);
radius_0=mean(sqrt(vcx(:,1).^2+vcx(:,2).^2+vcx(:,3).^2));

alpha=1;
rd = 1;
err_0=costfun(vc,center_0,radius_0);

for k=1:Ni;

  [center_new,radius_new]=lm1step(vc,center_0,radius_0,alpha);

  err_new=costfun(vc,center_new,radius_new);
  %     disp([k,err_0,err_new,center_new,radius_new]);

  if err_new<err_0;
    rd = abs(radius_0-radius_new)/radius_0;
    %         disp(rd)
    center_0=center_new;
    radius_0=radius_new;
    err_0=err_new;
    alpha=alpha/5;
  else
    alpha=alpha*5;
  end

  radius=radius_0;
  center=center_0;
  if rd<stopth, break, end % stop if

end

return;

%%
function err=costfun(vc,center,radius)
[nvc,ndum]=size(vc);

vcx=vc-repmat(center,nvc,1);

err=sqrt(sqrt(mean( (vcx(:,1).^2+vcx(:,2).^2+vcx(:,3).^2-radius^2).^2)));

return

%%
function  [center_new,radius_new]=lm1step(vc,center,radius,alpha)

[nvc,ndum]=size(vc);
vcx=vc-repmat(center,nvc,1);
f=vcx(:,1).^2+vcx(:,2).^2+vcx(:,3).^2-radius^2;

L=2*[vcx,repmat(radius,nvc,1)];

par_new=inv(L'*L+alpha*eye(4))*L'*f;

center_new=center+par_new(1:3)';
radius_new=radius+par_new(4);

return


