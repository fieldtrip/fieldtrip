function documentationreference(outdir);

funname = {
  'definetrial'
  'artifact_eog'
  'artifact_jump'
  'artifact_muscle'
  'rejectartifact'
  'preprocessing'
  'rejectvisual'
  'appenddata'
  'resampledata'
  'channelrepair'
  'recodeevent'
  'redefinetrial'
  'read_header'
  'read_data'
  'read_event'
  'read_mri'
  'read_spike'
  'write_data'
  'write_event'

  'timelockanalysis'
  'timelockgrandaverage'

  'freqanalysis'
  'freqanalysis_mtmfft'
  'freqanalysis_mtmwelch'
  'freqanalysis_mtmconvol'
  'freqanalysis_wltconvol'
  'freqanalysis_tfr'
  'freqbaseline'
  'freqgrandaverage'
  'freqdescriptives'

  'dipolefitting'
  'dipolesimulation'
  'sourceanalysis'
  'sourcegrandaverage'
  'sourcedescriptives'
  'sourceinterpolate'
  'prepare_localspheres'
  'prepare_singleshell'
  'prepare_bemmodel'
  'prepare_leadfield'
  'prepare_atlas'
  'volumelookup'
  'volumenormalise'
  'volumesegment'

  'timelockstatistics'
  'freqstatistics'
  'sourcestatistics'

  'spikedownsample'
  'spikesplitting'

  'rt_asaproxy'
  'rt_brainampproxy'
  'rt_ctfproxy'
  'rt_fileproxy'
  'rt_neuralynxproxy'
  'rt_plexonproxy'
  'rt_signalproxy'

  'prepare_layout'
  'layoutplot'
  'topoplot'
  'singleplotER'
  'singleplotTFR'
  'topoplotER'
  'topoplotTFR'
  'multiplotER'
  'multiplotTFR'
  'sourceplot'
  };

% create the desired output directory
if ~isdir(outdir)
mkdir(outdir);
end

for i=1:length(funname)
  filename = fullfile(outdir, [funname{i} '.txt']);
  str = help(funname{i});
  fid = fopen(filename, 'wt');
  fprintf(fid, '=====  %s =====\n\n', upper(funname{i}));
  fprintf(fid, 'Note that this reference documentation is identical to the help that is displayed in Matlab when you type "help %s".\n\n', funname{i});
  fprintf(fid, '<code>\n');  % required for docuwiki
  fprintf(fid, '%s', str);
  fprintf(fid, '</code>\n');  % required for docuwiki
  fclose(fid);
end


