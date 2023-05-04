function test_fieldtrip2fiff(mnedatadir)

%%
savedir = tempdir;

%%
% Use the mne-testing-data if path is specified
if nargin
  pw_dir = pwd;
  cd(mnedatadir);

%   datasets = {
%     fullfile(mnedatadir, 'CTF', 'catch-alp-good-f.ds')
%     fullfile(mnedatadir, 'CTF', 'somMDYO-18av.ds')
%     fullfile(mnedatadir, 'CTF', 'testdata_ctf.ds')
%     fullfile(mnedatadir, 'CTF', 'testdata_ctf_mc.ds')
%     fullfile(mnedatadir, 'CTF', 'testdata_ctf_pseudocontinuous.ds')
%     fullfile(mnedatadir, 'CTF', 'testdata_ctf_short.ds')
%     fullfile(mnedatadir, 'CTF', 'testdata_ctf_short_discontinuous.ds')
%     fullfile(mnedatadir, 'BTi', '4Dsim',   'c,rfDC')
%     fullfile(mnedatadir, 'BTi', 'erm_HFH', 'c,rfDC')
%     fullfile(mnedatadir, 'KIT', 'ArtificalSignalData_RICOH_1khz.con')
%     fullfile(mnedatadir, 'KIT', 'ArtificalSignalData_Yokogawa_1khz.con')
%     fullfile(mnedatadir, 'KIT', '010409_Motor_task_coregist-export_tiny_1s.con')
%     fullfile(mnedatadir, 'KIT', 'data_berlin.con')
%     fullfile(mnedatadir, 'MEG', 'sample', 'sample_audvis_trunc-ave.fif')
%     fullfile(mnedatadir, 'MEG', 'sample', 'sample_audvis-ave.fif')
%     };

  % the idea is here to test the status quo for reading and writing, hence
  % operate in first instance only on the fif-files
  d = dir('*/*.fif');
  d = cat(1, d, dir('*/*/*.fif'));
  d = cat(1, d, dir('*/*/*/*.fif'));
  
  datasets = cell(numel(d),1);
  for k = 1:numel(d)
    datasets{k} = fullfile(d(k).folder, d(k).name);
  end

  ok = false(numel(datasets),1);
  for k = 1:numel(datasets)
    cfg = [];
    cfg.dataset = datasets{k};
    try
      data = ft_preprocessing(cfg);
      ok(k,1) = true;
    catch
      ok(k,1) = false;
    end
  end


  cd(pw_dir);
end

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg');

% Use the 'Subject01.ds' dataset for basic testing of the writing and reading
cfg         = [];
cfg.dataset = fullfile(datadir, 'ctf151', 'Subject01.ds');
hdr         = ft_read_header(cfg.dataset);

% Define 2 epochs worth of data
cfg.trl     = [(0:1).*hdr.nSamples+1;(1:2).*hdr.nSamples;-300.*ones(1,2)]';
cfg.channel = 'MEG';
data        = ft_preprocessing(cfg);

% Benchmark trial for the original data
t1 = data.trial{1};

% Benchmark average for the original data
tlck = ft_timelockanalysis([], data);
a1 = tlck.avg;

% Save as epoched fif-files, either in double, or single precision format
fieldtrip2fiff(fullfile(savedir, 'data_epoched_double.fif'), data);
fieldtrip2fiff(fullfile(savedir, 'data_epoched_single.fif'), data, 'precision', 'single');

% Save as raw fif-files, either in double, or single precision format
cfgsel.trials = 1;
fieldtrip2fiff(fullfile(savedir, 'data_raw_double.fif'), ft_selectdata(cfgsel, data));
fieldtrip2fiff(fullfile(savedir, 'data_raw_single.fif'), ft_selectdata(cfgsel, data), 'precision', 'single');

% Save as evoked fif-files, either in double, or single precision format
fieldtrip2fiff(fullfile(savedir, 'data_evoked_double.fif'), ft_timelockanalysis([], data));
fieldtrip2fiff(fullfile(savedir, 'data_evoked_single.fif'), ft_timelockanalysis([], data), 'precision', 'single');

% Save as epoched/raw fif-files with complex data, either double or single precision
data.trial{1} = data.trial{1}+1i.*data.trial{2};
t1c  = data.trial{1};
tlck = ft_timelockanalysis([], data);
a1c  = tlck.avg;

fieldtrip2fiff(fullfile(savedir, 'data_epoched_complex_double.fif'), data);
fieldtrip2fiff(fullfile(savedir, 'data_epoched_complex_single.fif'), data, 'precision', 'single');
fieldtrip2fiff(fullfile(savedir, 'data_raw_complex_double.fif'), ft_selectdata(cfgsel, data));
fieldtrip2fiff(fullfile(savedir, 'data_raw_complex_single.fif'), ft_selectdata(cfgsel, data), 'precision', 'single');
fieldtrip2fiff(fullfile(savedir, 'data_evoked_complex_double.fif'), ft_timelockanalysis([], data));
fieldtrip2fiff(fullfile(savedir, 'data_evoked_complex_single.fif'), ft_timelockanalysis([], data), 'precision', 'single');

% Read in the data and compare against the benchmarks
cfg         = [];
cfg.dataset = fullfile(savedir, 'data_epoched_double.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, t1));
cfg.dataset = fullfile(savedir, 'data_epoched_single.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, single(t1)));
cfg.dataset = fullfile(savedir, 'data_epoched_complex_double.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, t1c));
cfg.dataset = fullfile(savedir, 'data_epoched_complex_single.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, single(t1c)));
cfg.dataset = fullfile(savedir, 'data_raw_double.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, t1));
cfg.dataset = fullfile(savedir, 'data_raw_single.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, single(t1)));
cfg.dataset = fullfile(savedir, 'data_raw_complex_double.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, t1c));
cfg.dataset = fullfile(savedir, 'data_raw_complex_single.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, single(t1c)));
cfg.dataset = fullfile(savedir, 'data_evoked_double.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, a1));
cfg.dataset = fullfile(savedir, 'data_evoked_single.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, single(a1)));
cfg.dataset = fullfile(savedir, 'data_evoked_complex_double.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, a1c));
cfg.dataset = fullfile(savedir, 'data_evoked_complex_single.fif');
datafif     = ft_preprocessing(cfg);
assert(isequal(datafif.trial{1}, single(a1c)));

%%
% This section tests a bunch of MEG dataset, of different systems
fname = {
  'bti148/c,rfhp0.1Hz'
  'bti248hcp/c,rfDC'
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

savename = cell(numel(fname),1);
for k = 1:numel(fname)
  cfg = [];
  cfg.dataset = fullfile(datadir, fname{k});
  if endsWith(cfg.dataset, 'fif')
    cfg.coilaccuracy = 0; %-> ensure the fif-reader to have the units right for the grad array, this uses the coil_def stuff
  end
  hdr = ft_read_header(cfg.dataset);
  if hdr.nTrials*hdr.nSamples>1000
    cfg.trl = [1 1000 0];
  end
  cfg.continuous = 'yes';
  data = ft_preprocessing(cfg);
  data.grad = ft_convert_units(data.grad, 'm');
  
  savename{k,1} = fullfile(savedir, sprintf('file%03d.fif',k));
  fieldtrip2fiff(savename{k}, data);
  save(strrep(savename{k},'fif','mat'),'data');
  clear data
end

for k = 1:numel(savename)
  skipgrad = false;

  cfg = [];
  cfg.dataset = savename{k};
  cfg.coilaccuracy = 0; % ensure fif readers to use mne2grad

  % the example data from the bti systems requires a custom coildef,
  % because of the representation of the diagonal reference gradiometers
  if contains(fname{k}, 'bti248/e')
    cfg.coildeffile = fullfile(datadir, 'bti248/coil_def_magnes_Glasgow.dat');
  elseif contains(fname{k}, 'bti248hcp')
    cfg.coildeffile = fullfile(datadir, 'bti248hcp/coil_def_magnes_StLouis.dat');
  elseif contains(fname{k}, 'bti248grad')
    cfg.coildeffile = fullfile(datadir, 'bti248grad/coil_def_magnes_Colorado.dat');
  elseif contains(fname{k}, 'A0132')
    % this one has a faulty bottom coil on one of the meggrads, will fail
    % for sure, not due to the functionality tested.
    skipgrad = true;
  end
  datafif = ft_preprocessing(cfg);
  load(strrep(savename{k},'fif','mat'));

  % some of the metadata are in single precision in the fif file
  data    = ft_struct2single(data);
  datafif = ft_struct2single(datafif);

  [ix,msg] = isalmostequal(rmfield(data,{'cfg' 'hdr' 'grad'}),rmfield(datafif,{'cfg' 'hdr' 'grad'}), 'reltol', 1e-4);
  M(k).msg = msg;

  if ~skipgrad
    % compare the grads
    grad    = data.grad;
    gradfif = datafif.grad;

    % the order of the channels might have been changed, as well as the
    % coils, as well as the polarity of the ori.
    [i1,i2] = match_str(grad.label, gradfif.label);
    fn = fieldnames(gradfif);
    for kk = 1:numel(fn)
      if size(gradfif.(fn{kk}), 1) == numel(i2)
        gradfif.(fn{kk}) = gradfif.(fn{kk})(i2,:);
      end
    end

    % reorder the coils
    i1 = zeros(1,0);
    i2 = zeros(1,0);
    for kk = 1:numel(grad.label)
      i1 = cat(2,i1,find(grad.tra(kk,:)~=0));
      i2 = cat(2,i2,find(gradfif.tra(kk,:)~=0));
    end

    for kk = 1:numel(fn)
      if size(gradfif.(fn{kk}), 1) == numel(i2)
        gradfif.(fn{kk}) = gradfif.(fn{kk})(i2,:);
      elseif size(gradfif.(fn{kk}), 2) == numel(i2)
        gradfif.(fn{kk}) = gradfif.(fn{kk})(:, i2);
      end
      if size(grad.(fn{kk}), 1) == numel(i1)
        grad.(fn{kk}) = grad.(fn{kk})(i1,:);
      elseif size(gradfif.(fn{kk}), 2) == numel(i1)
        grad.(fn{kk}) = grad.(fn{kk})(:, i1);
      end
    end

    sel = sum(gradfif.tra)<0;
    gradfif.tra = abs(gradfif.tra);
    gradfif.coilori(sel,:) = -gradfif.coilori(sel,:);
    
    if endsWith(fname{k}, 'fif')
      % this does not happen too often: but for the native fif-files
      sel = sum(grad.tra)<0;
      grad.tra = abs(grad.tra);
      grad.coilori(sel,:) = -grad.coilori(sel,:);
    end

    assert(all(sum(grad.coilori.*gradfif.coilori,2)>0.999), 'coil orientation different');
    assert(all(sqrt(sum((grad.coilpos-gradfif.coilpos).^2,2))<1e-3), 'coil position different');
  end
end
