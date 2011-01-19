function [Y,LY] = hh_downsampling(gridlocs,localmax,voxel_sizes,min_dis)

%% Init parameters
X = double(gridlocs);
LX = localmax;
Y = [];
LY = [];
MINDIS = min_dis * min_dis;

%% Downsampling the given set of nodes
while ~isempty(X)
    Y = [Y;X(1,:)];
    LY = [LY,LX(1)];
    dx = (X(:,1) - X(1,1)) * voxel_sizes(1);
    dy = (X(:,2) - X(1,2)) * voxel_sizes(2);
    dz = (X(:,3) - X(1,3)) * voxel_sizes(3);
    dis = dx.^2 + dy.^2 + dz.^2;
    idx = dis > MINDIS;
    X = X(idx,:);
    LX = LX(idx);
end
