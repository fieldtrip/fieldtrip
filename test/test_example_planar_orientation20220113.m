function test_example_planar_orientation

% MEM 4gb
% WALLTIME 00:10:00

%
%% Find the orientation of planar gradiometers
%
% The Neuromag/Elekta/Megin Vectorview and Triux scanners have a combination of magnetometers and two planar gradiometers at each sensor location.
%
% This is further detailed in the Neuromag/Elekta/MEGIN System Hardware Users' Guide and on [this page](http://imaging.mrc-cbu.cam.ac.uk/meg/VectorviewDescription#Magsgrads).
%
% The following code demonstrates how you can determine the orientation of the two planar channels at each location.
%

filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/workshop/natmeg/oddball1_mc_downsampled.fif');
grad = ft_read_sens(filename, 'senstype', 'meg');

sel = find(strcmp(grad.chantype, 'megplanar'));

for i=1:numel(sel)
j = sel(i);
clear coilindex*
coilindex1 = find(grad.tra(j,:)>0)
coilindex2 = find(grad.tra(j,:)<0)
pos(i,:) = (grad.coilpos(coilindex2,:) + grad.coilpos(coilindex1,:))/2;
ori(i,:) = grad.coilpos(coilindex2,:) - grad.coilpos(coilindex1,:);
ori(i,:) = ori(i,:)/norm(ori(i,:));

end

figure
quiver3(pos(:,1), pos(:,2), pos(:,3), ori(:,1), ori(:,2), ori(:,3))
axis equal
axis vis3d
