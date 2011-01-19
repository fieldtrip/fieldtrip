clear;

%% Init the test data provided by MATLAB
[x,y,z,v] = flow(10); 
[xi,yi,zi] = meshgrid(.1:.25:10,-3:.25:3,-3:.25:3);
vi = interp3(x,y,z,v,xi,yi,zi); % vi is 25-by-40-by-25
slice(xi,yi,zi,vi,[6 9.5],2,[-2 .2]), shading flat;

%% Display the test data
scatter3(x(:),y(:),z(:),abs(v(:)))

%% Test with the tomographic map