%% Parameters
addpath /home/hungptit/working/fdm/externals/bioelectromagnetism;

%% Obtain the scalp
unix('rm scalp.h5');
unix('../../bin/img_get_scalp --imgfile smri-seghead --rgnfile scalp.h5');
scalp = hh_region_read('scalp.h5');
plot3(scalp.gridlocs(:,1),scalp.gridlocs(:,2),scalp.gridlocs(:,3),'*');

%% Read in electrode locations
elecs = importdata('electrode.txt');
nas = elecs.data(1,:);
lpa = elecs.data(2,:);
rpa = elecs.data(3,:);
P = [nas;rpa;lpa];
NAS = [87,209,141];
RPA = [164,98,120];
LPA = [13,99,120];
Q = [NAS;RPA;LPA];

%% Compute the rotation matrix
T = fiducial_coregister(P,Q);
P2Q = [P,ones(size(P,1),1)] * T;
X = elecs.data;

%% Fit electrode to the scalp
electrodes = hh_elecs2scalp(X,scalp,P,Q);

%% Save the result to the output file
electrodes.info = 'Fitted electrodes';
hh_region_write('electrodes.h5',electrodes);
unix(['../../bin/rgn_img --img smri-seghead --rgn ' ...
      'electrodes.h5 --rgnimg electrodes']);

%% Obtain the reciprocity data
unix(['elecsfwd --img smri-seghead --electrodes electrodes.h5 --data ' ...
      'data.h5 --contable ../../data/conductivity/contable.csv --TOL 1e-8']);


%% Extract the reciprocity data for the brain region only

%% Generate the simulated sources

%% Obtain the output map of LCMV

%% Test the raw results
% $$$ NELLI = size(gridlocs,1);
% $$$ electrodes.node_sizes = scalp.node_sizes;
% $$$ electrodes.voxel_sizes = scalp.voxel_sizes;
% $$$ electrodes.info = 'Fitted electrodes';
% $$$ electrodes.gridlocs = int32(round(gridlocs));
% $$$ electrodes.locations = gridlocs;
% $$$ electrodes.values = zeros(NELEC,1);
% $$$ electrodes.status = ones(NELEC,1)
% $$$ electrodes.info = 'Fitted electrodes';
% $$$ hh_region_write('electrodes.h5',electrodes);

%% Display the data using 3D plot
% $$$ plot3(scalp.gridlocs(:,1),scalp.gridlocs(:,2),scalp.gridlocs(:,3),'r*');
% $$$ hold on;
% $$$ plot3(gridlocs(:,1),gridlocs(:,2),gridlocs(:,3),'*');
% $$$ hold off;
