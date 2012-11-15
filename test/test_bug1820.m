function test_bug1820

% TEST test_bug1820
% TEST ft_headmodel_simbio ft_prepare_vol_sens ft_compute_leadfield

% See http://bugzilla.fcdonders.nl/show_bug.cgi?id=1820

cd /home/common/matlab/fieldtrip/data/test/bug1820

% load the mesh
load test_forward

% first calculate stiffness matrix, will be stored in test_stiff.stiff
test_stiff = ft_headmodel_simbio(test,'conductivity',[0.33,0.0042,0.33]);

% now compute transfer matrix, will be stored in test_stiff.transfer
test_transfer = ft_prepare_vol_sens(test_stiff,sens);

% finally compute the leadfield
lf = ft_compute_leadfield(pos,sens,test_transfer);

% now it would be good to compute leadfields with an analytic concentric
% sphere solution and compare them
