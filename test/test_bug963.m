function test_bug963

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_read_header ft_read_sens ft_datatype_sens bti2grad itab2grad netmeg2grad ctf2grad mne2grad yokogawa2grad fif2grad mne2grad.old yokogawa2grad_new ft_compute_leadfield ft_prepare_vol_sens

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/bug963');
rawdataprefix = dccnpath('/home/common/matlab/fieldtrip/data/test');

dataset = {
  'original/meg/bti148/c,rfhp0.1Hz'
  'original/meg/bti248/e,rfDC,F,a'
  'original/meg/bti248grad/e,rfhp1.0Hz,COH'
  'original/meg/ctf151/Subject01.ds'
  'original/meg/ctf275/A0132_Aud-Obj-Recognition_20051115_02.ds'
  'original/meg/ctf275/TacStimRegressConfound.ds'
  'original/meg/ctf64/Wat123r1raw.ds'
  'original/meg/itab153/srgcst85_0105.raw.mhd'
  'original/meg/itab28/gibb0101.raw.mhd'
  'original/meg/itab28_old/gibb0101.raw.mhd'
  'original/meg/neuromag122/jg_single_01raw.fif'
  'original/meg/neuromag306/raw.fif'
  'original/meg/yokogawa160/Continuous1.con'
  'original/meg/yokogawa440/S1_MEG_Epoch.raw'
  'original/meg/yokogawa64/2011_01_28_0354_ME053_AEF.con'
  };

%cd(datadir)

for i=1:length(dataset)
  
  if i==7
    % for this one the autodetection is not fine-grained enough
    headerformat = 'ctf_old';
  else
    headerformat = [];
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 1A: generate the reference data as a set of *.mat files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  outputfile = fullfile(datadir, sprintf('dataset%02d.mat', i));
  if ~exist(outputfile)
    filename = fullfile(rawdataprefix, dataset{i});
    disp(filename)
    hdr  = ft_read_header(filename, 'headerformat', headerformat);
    grad = ft_read_sens(filename, 'fileformat', headerformat, 'senstype', 'meg');
    warning('writing reference solution to %s', outputfile);
    save(outputfile, 'hdr', 'grad');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART 1B: read the datasets and compare them to the reference data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  outputfile = fullfile(datadir, sprintf('dataset%02d.mat', i));
  reference = load(outputfile);
  
  filename = fullfile(rawdataprefix, dataset{i});
  disp(filename)
  hdr  = ft_read_header(filename, 'headerformat', headerformat);
  grad = ft_read_sens(filename, 'fileformat', headerformat, 'senstype', 'meg');
  
  % remove the grad.balance field if the current balancing is none
  % as that it is not interesting to compare
  if isfield(grad, 'balance') && strcmp(grad.balance.current, 'none')
    grad = rmfield(grad, 'balance');
  end
  if isfield(hdr.grad, 'balance') && strcmp(hdr.grad.balance.current, 'none')
    hdr.grad = rmfield(hdr.grad, 'balance');
  end
  if isfield(reference.grad, 'balance') && strcmp(reference.grad.balance.current, 'none')
    reference.grad = rmfield(reference.grad, 'balance');
  end
  if isfield(reference.hdr.grad, 'balance') && strcmp(reference.hdr.grad.balance.current, 'none')
    reference.hdr.grad = rmfield(reference.hdr.grad, 'balance');
  end
  
  % remove coordsys field as these were not yet present in reference files
%   if isfield(grad, 'coordsys')
%     grad = rmfield(grad, 'coordsys');
%   end
%   if isfield(hdr.grad, 'coordsys')
%     hdr.grad = rmfield(hdr.grad, 'coordsys');
%   end
%   if isfield(grad, 'labelold')
%     grad = rmfield(grad, 'labelold');
%   end
%   if isfield(hdr.grad, 'labelold')
%     hdr.grad = rmfield(hdr.grad, 'labelold');
%   end
%   
  assert(isalmostequal(hdr.grad,           grad, 'reltol',eps*1e6), sprintf('failed for %s', filename));
  assert(isalmostequal(reference.grad,     grad, 'reltol',eps*1e6), sprintf('failed for %s', filename));
  assert(isalmostequal(reference.hdr.grad, grad, 'reltol',eps*1e6), sprintf('failed for %s', filename));
  
  allhdr{i}  = hdr;
  allgrad{i} = grad;
  
end % for all datasets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 2:
% also do a quick sanity check on the tutorial data
% this checks that all three CTF implementations still work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
hdr1 = ft_read_header(filename, 'headerformat', 'ctf_ds');
hdr2 = ft_read_header(filename, 'headerformat', 'read_ctf_res4');
%hdr3 = ft_read_header(filename, 'headerformat', 'ctf_read_res4');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 3:
% it is important that all missing information can be added automatically
% since people might have old grad structures in *.mat files on disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(dataset)
  
  grad1  = allgrad{i};
  
  grad2 = grad1;
  try
    grad2 = rmfield(grad2, 'type');
  end
  grad2.type = ft_senstype(grad2);
  assert(~any(strcmp(grad2.type, {'unknown', 'meg'})));
  
  grad3 = grad1;
  try
    grad3 = rmfield(grad3, 'chantype');
  end
  grad3.chantype = ft_chantype(grad3);
  assert(mean(strcmp('unknown', grad3.chantype))<0.3);
  
  grad4 = grad1;
  try
    grad4 = rmfield(grad4, 'chanunit');
  end
  grad4.chanunit = ft_chanunit(grad4);
  assert(mean(strcmp('unknown', grad4.chanunit))<=mean(strcmp('unknown', grad1.chanunit))); % the number of channels with unknown units should not increase
  
  grad5 = grad1;
  try
    grad5 = rmfield(grad5, 'type');
  end
  try
    grad5 = rmfield(grad5, 'chantype');
  end
  try
    grad5 = rmfield(grad5, 'chanunit');
  end
  grad5 = ft_datatype_sens(grad5);
  assert(~any(strcmp(grad5.type, {'unknown', 'meg'})));
  assert(mean(strcmp('unknown', grad5.chanunit))<=mean(strcmp('unknown', grad1.chanunit))); % the number of channels with unknown units should not increase
  assert(mean(strcmp('unknown', grad5.chantype))<=mean(strcmp('unknown', grad1.chantype))); % the number of channels with unknown type should not increase
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 4:
% try out the unit conversion on a hand-crafted gradiometer array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sens          = [];
sens.coilpos  = [0.01 0 0.01; -0.01 0 0.01];
sens.coilori  = [0 0 1; 0 0 1];
sens.tra      = [1 -1] / 0.02; % divide by the baseline
sens.label    = {'M1'};
sens.chantype = {'megplanar'};
sens.chanunit = {'T/m'};
sens.unit     = 'm';
sens.type     = 'meg';

if false
  sens          = [];
  sens.coilpos  = [0 0 0.01];
  sens.coilori  = [0 0 1];
  sens.chanpos  = [0 0 0.01];
  sens.chanori  = [0 0 1];
  sens.tra      = 1;
  sens.label    = {'M1'};
  sens.chantype = {'megmag'};
  sens.chanunit = {'T'};
  sens.unit     = 'm';
  sens.type     = 'meg';
end

vol      = [];
vol.type = 'infinite';
vol.unit = 'm';

pos = [0 0 -0.01];

pos_m  = pos;
pos_cm = pos*100;
pos_mm = pos*1000;

vol_m  = ft_convert_units(vol, 'm');
vol_cm = ft_convert_units(vol, 'cm');
vol_mm = ft_convert_units(vol, 'mm');

sens_m  = ft_convert_units(sens, 'm');
sens_cm = ft_convert_units(sens, 'cm');
sens_mm = ft_convert_units(sens, 'mm');

if strcmp(sens.chantype{1}, 'megplanar')
  % this should be scaled with the geometrical units
  assert(~isequal(sens_m.tra, sens_cm.tra));
  assert(~isequal(sens_m.tra, sens_cm.tra));
end

