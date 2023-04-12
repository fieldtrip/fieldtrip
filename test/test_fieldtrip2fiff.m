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
  cfg = [];
  cfg.dataset = savename{k};
  cfg.coilaccuracy = 0; % ensure fif readers to use mne2grad
  datafif = ft_preprocessing(cfg);
  load(strrep(savename{k},'fif','mat'));

  % some of the metadata are in single precision in the fif file
  data    = ft_struct2single(data);
  datafif = ft_struct2single(datafif);

  [ix,msg] = isalmostequal(rmfield(data,{'cfg' 'hdr' 'grad'}),rmfield(datafif,{'cfg' 'hdr' 'grad'}), 'reltol', 1e-4);
  M(k).msg = msg;

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
  
  c1 = grad.coilpos;
  c2 = gradfif.coilpos;
  n  = size(c1,1);
  D  = squareform(pdist([c1;c2]));
  D  = D(1:n,n+(1:n));
  
  i2 = zeros(size(D,1),1);
  for kk = 1:numel(i2)
    [m, i2(kk)] = min(D(kk,:));
  end
  
  for kk = 1:numel(fn)
    if size(gradfif.(fn{kk}), 1) == numel(i2)
      gradfif.(fn{kk}) = gradfif.(fn{kk})(i2,:);
    elseif size(gradfif.(fn{kk}), 2) == numel(i2)
      gradfif.(fn{kk}) = gradfif.(fn{kk})(:, i2);
    end
  end
  
  sel = sum(gradfif.tra)<0;
  gradfif.tra = abs(gradfif.tra);
  gradfif.coilori(sel,:) = -gradfif.coilori(sel,:);

end
