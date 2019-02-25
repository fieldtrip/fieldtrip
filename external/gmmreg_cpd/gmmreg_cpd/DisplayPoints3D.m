function [axis_limits] = DisplayPoints3D(Model, Scene, sampling, axis_limits)
%%=====================================================================
%% $RCSfile: DisplayPoints3D.m,v $
%% $Date: 2008/12/05 00:09:32 $
%% $Revision: 1.1 $
%%=====================================================================

set(gca,'FontSize',16,'FontName','Times','FontWeight','bold');

plot3(Model(:,1),Model(:,2),Model(:,3),'r+', 'MarkerSize', 1.5, 'LineWidth',1);
hold on;
plot3(Scene(:,1),Scene(:,2),Scene(:,3),'bo', 'MarkerSize', 1.5, 'LineWidth',1);
axis equal;

if (nargin<3)
%    axis_limits = determine_border(Model, Scene);
    sampling = 0;
end

m = size(Model,1);
if (sampling>0)
    for i=1:sampling:m
        text(Model(i,1), Model(i,2), Model(i,3), [' \leftarrow',sprintf('%d',i)]);
    end
end

m = size(Scene,1);
if (sampling>0)
    for i=1:sampling:m
        text(Scene(i,1), Scene(i,2), Scene(i,3), [' \leftarrow',sprintf('%d',i)]);
    end
end

if (nargin<4)
    axis_limits = determine_border(Model, Scene);
end

xlim(axis_limits(1,:));
ylim(axis_limits(2,:));   
zlim(axis_limits(3,:));   
axis off
%pbaspect([1,1,1]);