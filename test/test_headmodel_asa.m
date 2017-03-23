function test_headmodel_asa

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_headmodel_asa ft_prepare_vol_sens ft_compute_leadfield

% the standard BEM model that is supplied on the FieldTrip ftp server 
% and included in EEGLAB was originally constructed in ASA version 2

% see ftp://ftp.fcdonders.nl/pub/fieldtrip/standard_BEM.zip
% although that zip file does not contain the original file ASA formats

hdmfile = dccnpath('/home/common/matlab/fieldtrip/template/headmodel/standard_bem.mat');
elcfile = dccnpath('/home/common/matlab/fieldtrip/template/electrode/standard_1020.elc');

vol = ft_headmodel_asa(hdmfile);
figure
ft_plot_vol(vol, 'facecolor', 'none')

if true
  % read some standard electrode positions
  sens = ft_read_sens(elcfile);
else
  % create a random set of electrodes
  pnt = randn(200,3);
  for i=1:200
    pnt(i,:) = pnt(i,:) ./ norm(pnt(i,:))
  end
  sel = find(pnt(:,3)>0);
  sens = [];
  sens.elecpos = pnt(sel,:) * 100;
  sens.chanpos = pnt(sel,:) * 100;
  for i=1:length(sel)
    sens.label{i} = sprintf('chan%03d', i);
  end
end

% the electrodes will be projected on the surface
orig = sens;

% project the electrodes on the volume conduction model
[vol, sens] = ft_prepare_vol_sens(vol, sens);


% compare the original and projected electrode positions
figure
ft_plot_sens(sens);
ft_plot_sens(orig, 'style', 'ro');

% conpute an example leadfield
lf = ft_compute_leadfield([0 0 50], sens, vol);

figure;
subplot(2,2,1); ft_plot_topo3d(sens.chanpos, lf(:,1))
subplot(2,2,2); ft_plot_topo3d(sens.chanpos, lf(:,2))
subplot(2,2,3); ft_plot_topo3d(sens.chanpos, lf(:,3))



