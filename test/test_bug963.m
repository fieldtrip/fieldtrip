function test_bug963

% TEST test_bug963
% TEST ft_read_header ft_read_sens ft_datatype_sens bti2grad itab2grad netmeg2grad ctf2grad mne2grad yokogawa2grad fif2grad mne2grad.old yokogawa2grad_new

datadir = '/home/common/matlab/fieldtrip/data/test/bug963';
rawdataprefix = '/home/common/matlab/fieldtrip/data/test';

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

cd(datadir)

for i=1:length(dataset)
  
  if i==7
    % for this one the autodetection is not fine-grained enough
    headerformat = 'ctf_old';
  else
    headerformat = [];
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART ONE: generate the reference data as a set of *.mat files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  outputfile = sprintf('dataset%02d.mat', i);
  if ~exist(outputfile)
    filename = fullfile(rawdataprefix, dataset{i});
    disp(filename)
    hdr  = ft_read_header(filename, 'headerformat', headerformat);
    grad = ft_read_sens(filename, 'fileformat', headerformat);
    save(outputfile, 'hdr', 'grad');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PART TWO: read the datasets and compare them to the reference data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  outputfile = sprintf('dataset%02d.mat', i);
  reference = load(outputfile);
  
  filename = fullfile(rawdataprefix, dataset{i});
  disp(filename)
  hdr  = ft_read_header(filename, 'headerformat', headerformat);
  grad = ft_read_sens(filename, 'fileformat', headerformat);
  
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
    
  assert(isequal(hdr.grad,           grad), sprintf('failed for %s', filename));
  assert(isequal(reference.grad,     grad), sprintf('failed for %s', filename));
  assert(isequal(reference.hdr.grad, grad), sprintf('failed for %s', filename));
  
  allhdr{i}  = hdr;
  allgrad{i} = grad;
  
end % for all datasets


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% also do a quick sanity check on the tutorial data
% this checks that all three CTF implementations still work

filename = '/home/common/matlab/fieldtrip/data/Subject01.ds';
hdr1 = ft_read_header(filename, 'headerformat', 'ctf_ds');
hdr2 = ft_read_header(filename, 'headerformat', 'read_ctf_res4');
hdr3 = ft_read_header(filename, 'headerformat', 'ctf_read_res4');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% it is important that all missing information can be added automatically
% since people might have old grad structures in *.mat files on disk

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


