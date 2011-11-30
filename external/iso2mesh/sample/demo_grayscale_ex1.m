% this example demonstrates how to mesh a gray-scale volume directly.
% the data set is a Jacobian matrix obtained by multiplying two
% Green's functions
% see Appendix in Qianqian Fang, Computational methods for microwave medical imaging. Ph.D.
% Dissertation, Dartmouth College, Hanover, NH, U.S.A., Dec 2004.

[xi,yi,zi]=meshgrid(1:40,1:40,1:80);
r1=[xi(:)-20,yi(:)-20,zi(:)-20];
r2=[xi(:)-20,yi(:)-20,zi(:)-60];
r1=reshape(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2,size(xi));
r2=reshape(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2,size(xi));

k=10;

% generate Green's functions

g1=exp(sqrt(-1)*k*r1)./(4*pi*r1);
g2=exp(sqrt(-1)*k*r2)./(4*pi*r2);


%get Jacobian (sensitivity)
g12=g1.*g2;

% create a mesh based on the level-set at log(J)=-7.5

[node,elem,face]=v2m(log10(abs(g12))+10,2.5,2,10,'cgalsurf');

plotmesh(node,face);
