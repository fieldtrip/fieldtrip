function test_fieldtrip2fiff

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg');

fname = {
  'bti148/c,rfhp0.1Hz'
  'bti248/e,rfDC'
  'bti248grad/e,rfhp1.0Hz,COH'
  'ctf275/A0132_Aud-Obj-Recognition_20051115_02.ds'
  'ctf275/TacStimRegressConfound.ds'
  'ctf151/Subject01.ds'
  'neuromag122/jg_single_01raw.fif'
  'neuromag306/raw.fif'
  'neuromag306/run_01_raw.fif'
  'neuromag306/sub-15_ses-meg_task-facerecognition_run-01_meg.fif'
  };

savedir = tempdir;
savename = cell(numel(fname),1);
for k = 1:numel(fname)
  cfg = [];
  cfg.dataset = fullfile(datadir, fname{k});
  cfg.coilaccuracy = 0; %-> ensure the fif-reader to have the units right for the grad array
  hdr = ft_read_header(cfg.dataset);
  if hdr.nTrials*hdr.nSamples>1000
    cfg.trl = [1 1000 0];
  end
  cfg.continuous = 'yes';
  data = ft_preprocessing(cfg);
  
  savename{k,1} = fullfile(savedir, sprintf('file%03d.fif',k));
  fieldtrip2fiff(savename{k}, data);
  save(strrep(savename{k},'fif','mat'),'data');
  clear data
end

for k = 1:numel(savename)
  cfg = [];
  cfg.dataset = savename{k};
  cfg.coilaccuracy = 0; % ensure fif readers to use mne2grad
  datafif = ft_preprocessing(cfg);
  load(strrep(savename{k},'fif','mat'));

  % some of the metadata are in single precision in the fif file
  data    = ft_struct2single(data);
  datafif = ft_struct2single(datafif);

  [ix,msg] = isalmostequal(rmfield(data,{'cfg' 'hdr'}),rmfield(datafif,{'cfg' 'hdr'}));
end
