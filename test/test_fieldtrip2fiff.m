function test_fieldtrip2fiff

savedir = tempdir;

% test a variety of MEG files
datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg');

cfg = [];
cfg.dataset = fullfile(datadir, 'ctf151', 'Subject01.ds');
hdr = ft_read_header(cfg.dataset);

cfg.trl = [(0:4).*hdr.nSamples+1;(1:5).*hdr.nSamples;-300.*ones(1,5)]';
data = ft_preprocessing(cfg);
fieldtrip2fiff(fullfile(savedir, 'data_epoched.fif'), data);

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
