function test_bug3172

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_plot_sens

close all

%% EEG electrodes

sens = ft_read_sens('easycap-M10.txt');
figure; ft_plot_sens(sens);
figure; ft_plot_sens(sens, 'coilshape', 'point');
if false
  % these two should error
  figure; ft_plot_sens(sens, 'coilshape', 'circle');
  figure; ft_plot_sens(sens, 'coilshape', 'square');
end

% giving it orientation allows the electrodes to be plotted as circles
sens.chanori = sens.chanpos;
sens.elecori = sens.elecpos;
for i=1:61
  sens.chanori(i,:) = sens.chanori(i,:) ./ norm(sens.chanori(i,:));
  sens.elecori(i,:) = sens.elecori(i,:) ./ norm(sens.elecori(i,:));
end

figure; ft_plot_sens(sens, 'coilshape', 'circle');

if false
  % this one should error
  figure; ft_plot_sens(sens, 'coilshape', 'square');
end


%% CTF151

sens = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf151/Subject01.ds'));
figure; ft_plot_sens(sens);
figure; ft_plot_sens(sens, 'coilshape', 'point', 'coil', true);
figure; ft_plot_sens(sens, 'coilshape', 'point', 'coil', false);
figure; ft_plot_sens(sens, 'coilshape', 'circle', 'coil', true);
figure; ft_plot_sens(sens, 'coilshape', 'circle', 'coil', false);
if false
  % this one should error
  figure; ft_plot_sens(sens, 'coilshape', 'square', 'coil', true);
  figure; ft_plot_sens(sens, 'coilshape', 'square', 'coil', false);
end

%% NEUROMAG306

sens = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/raw.fif'), 'senstype', 'meg');
figure; ft_plot_sens(sens);
figure; ft_plot_sens(sens, 'chantype', 'megmag');
figure; ft_plot_sens(sens, 'chantype', 'megplanar');
figure; ft_plot_sens(sens, 'coilshape', 'point', 'coil', true);
figure; ft_plot_sens(sens, 'coilshape', 'point', 'coil', false);
figure; ft_plot_sens(sens, 'coilshape', 'circle', 'coil', true, 'coilsize', 0.2); % small circles
figure; ft_plot_sens(sens, 'coilshape', 'circle', 'coil', false);
if false
  % this one should error
  figure; ft_plot_sens(sens, 'coilshape', 'square', 'coil', true);
end
figure; ft_plot_sens(sens, 'coilshape', 'square', 'coil', false);

figure; ft_plot_sens(sens, 'coilshape', 'square', 'coil', false, 'chantype', 'megmag'); % should be empty figure
figure; ft_plot_sens(sens, 'coilshape', 'square', 'coil', false, 'chantype', 'megplanar'); % should not be empty figure

%% NEUROMAG122

sens = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag122/jg_single_01raw.fif'));

figure; ft_plot_sens(sens);
figure; ft_plot_sens(sens, 'coilshape', 'point', 'coil', true);
figure; ft_plot_sens(sens, 'coilshape', 'point', 'coil', false);
figure; ft_plot_sens(sens, 'coilshape', 'circle', 'coil', true, 'coilsize', 0.2); % small circles
figure; ft_plot_sens(sens, 'coilshape', 'circle', 'coil', false);
if false
  % this one should error
  figure; ft_plot_sens(sens, 'coilshape', 'square', 'coil', true);
end
figure; ft_plot_sens(sens, 'coilshape', 'square', 'coil', false);

