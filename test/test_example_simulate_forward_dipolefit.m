function test_example_simulate_forward_dipolefit

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_dipolesimulation ft_timelockanalysis ft_dipolefitting

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% create a set of electrodes, randomly placed on the sphere
elec = [];
elec.pnt = randn(128,3);
dum = sqrt(sum(elec.pnt.^2,2));
elec.pnt = elec.pnt ./ [dum dum dum];  % scale them to a unit sphere
for i=1:128
   elec.label{i} = sprintf('%03d', i);
end

% Note: I suggest that you do not use these electrodes, but instead that 
% you get more meaningful electrode locations. See for example  
% http://robertoostenveld.ruhosting.nl/index.php/electrode/

% create a concentric 3-sphere volume conductor, the radius is the same as for the electrodes
vol = [];
vol.r = [0.88 0.92 1.00]; % radii of spheres
vol.cond = [1 1/80 1];       % conductivity
vol.o = [0 0 0];          % center of sphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% create a dipole simulation with one dipole and a 10Hz sine wave
cfg      = [];
cfg.vol  = vol;             % see above
cfg.elec = elec;            % see above
cfg.dip.pos = [0 0.5 0.3];
cfg.dip.mom = [1 0 0]';     % note, it should be transposed
cfg.dip.frequency = 10;
cfg.ntrials = 10;
cfg.triallength = 1;        % seconds
cfg.fsample = 250;          % Hz
raw1 = ft_dipolesimulation(cfg);
avg1 = ft_timelockanalysis([], raw1);
plot(avg1.time, avg1.avg);  % plot the timecourse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% create a dipole simulation with one dipole and a custom timecourse
cfg      = [];
cfg.vol  = vol;               % see above
cfg.elec = elec;              % see above
cfg.dip.pos = [0 0.5 0.3];
cfg.dip.mom = [1 0 0]';       % note, it should be transposed
cfg.fsample = 250;            % Hz
time = (1:250)/250;           % manually create a time axis
signal = sin(10*time*2*pi);   % manually create a signal
cfg.dip.signal = {signal, signal, signal};  % three trials
raw2 = ft_dipolesimulation(cfg);
avg2 = ft_timelockanalysis([], raw2);
plot(avg2.time, avg2.avg);    % plot the timecourse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% create a dipole simulation with two dipoles and a custom timecourse
cfg      = [];
cfg.vol  = vol;      % see above
cfg.elec = elec;     % see above
cfg.dip.pos = [
   0  0.5 0.3        % dipole 1
   0 -0.5 0.3        % dipole 2
   ];
cfg.dip.mom = ...       % the vector represents [qx1 qy1 qz1 qx2 qy2 qz2]
  [ 1 0 0 0 0 0 ]' + ...% this is how signal1 contributes to the 6 dipole components
  [ 0 0 0 1 0 0 ]';     % this is how signal2 contributes to the 6 dipole components
                        % note, it should be transposed
time = (1:250)/250;
signal1 = sin(10*time*2*pi);
signal2 = cos(15*time*2*pi);
cfg.dip.signal = {[signal1; signal2]}; % one trial only
cfg.fsample = 250;                     % Hz
raw3 = ft_dipolesimulation(cfg);
avg3 = ft_timelockanalysis([], raw3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% do a dipole fit of the first simulated dataset
cfg      = [];
cfg.vol  = vol;         % see above
cfg.elec = elec;        % see above
cfg.dip.pos = [0 0 0];  % initial search position
cfg.gridsearch = 'no';
dip1 = ft_dipolefitting(cfg, avg1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% or alternatively start with an exhaustive search on a coarse grid
cfg = [];
cfg.vol  = vol;         % see above
cfg.elec = elec;        % see above
cfg.gridsearch = 'yes';
cfg.xgrid = linspace(-1,1,5);
cfg.ygrid = linspace(-1,1,5);
cfg.zgrid = linspace(-1,1,5);
dip2 = ft_dipolefitting(cfg, avg1);
