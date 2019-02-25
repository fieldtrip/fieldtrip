function [axis_limits] = determine_border(Model, Scene)
%%=====================================================================
%% Date:      $Date: 2008/12/05 00:09:32 $
%% Version:   $Revision: 1.1 $
%%=====================================================================

dim = size(Scene,2);
axis_limits = zeros(dim,2);
for i=1:dim
    min_i = min([Scene(:,i);Model(:,i)]);
    max_i = max([Scene(:,i);Model(:,i)]);
    margin_i = (max_i-min_i)*0.05;
    axis_limits(i,:) = [min_i - margin_i max_i+margin_i];
end

