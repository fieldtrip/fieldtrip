function failed_bug2474

% MEM 1000mb
% WALLTIME 00:10:00

% TEST test_bug2474
% TEST ft_compute_leadfield ft_prepare_vol_sens


load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2474/bug.mat'));

try
  % this initially gives an error
  lf = ft_compute_leadfield(pos, sens, vol, 'reducerank', 2, 'dipoleunit', 'nA*m', 'chanunit', chanunits);
catch
  disp('the expected error happened');
end

% this should not give an error
lf = ft_compute_leadfield(pos, ft_struct2double(sens), ft_struct2double(vol), 'reducerank', 2, 'dipoleunit', 'nA*m', 'chanunit', chanunits);


% from 2014 onwards, these have to be in double precision
sens = ft_datatype_sens(sens, 'version', 'upcoming');
vol  = ft_datatype_headmodel(vol);

% this should now not be a problem any more
lf = ft_compute_leadfield(pos, sens, vol, 'reducerank', 2, 'dipoleunit', 'nA*m', 'chanunit', chanunits);
