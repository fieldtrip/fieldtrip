function test_bug168

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_realtime_topography
 
[ftver, ftpath] = ft_version;
cd(ftpath);

% ensure that these compat folders are presaent on the path
% although they should not be used
ft_hastoolbox('COMPAT', 2);
ft_hastoolbox('UTILITIES/COMPAT', 2);
ft_hastoolbox('FILEIO/COMPAT', 2);
ft_hastoolbox('PREPROC/COMPAT', 2);
ft_hastoolbox('FORWARD/COMPAT', 2);
ft_hastoolbox('PLOTTING/COMPAT', 2);
    
% this list was compiled on 1 May 2012 from the fieldtrip/realtime directory
% functions that will be added future don't have to be tested persee
inlist = {
  'acquisition/matlab/ft_realtime_asaproxy.m'
  'acquisition/matlab/ft_realtime_brainampproxy.m'
  'acquisition/matlab/ft_realtime_ctfproxy.m'
  'acquisition/matlab/ft_realtime_dicomproxy.m'
  'acquisition/matlab/ft_realtime_fileproxy.m'
  'acquisition/matlab/ft_realtime_fmriproxy.m'
  'acquisition/matlab/ft_realtime_neuralynxproxy.m'
  'acquisition/matlab/ft_realtime_pooraudioproxy.m'
  'acquisition/matlab/ft_realtime_signalproxy.m'
  'acquisition/matlab/private/encode_nifti1.m'
  'acquisition/matlab/private/tcpread.m'
  'acquisition/siemens/compile.m'
  'acquisition/siemens/matlab/replay_dicoms.m'
  'buffer/matlab/buffer.m'
  'buffer/matlab/compile.m'
  'example/ft_realtime_average.m'
  'example/ft_realtime_benchmark.m'
  'example/ft_realtime_classification.m'
  'example/ft_realtime_downsample.m'
  'example/ft_realtime_fmriviewer.m'
  'example/ft_realtime_heartratemonitor.m'
  'example/ft_realtime_packettimer.m'
  'example/ft_realtime_powerestimate.m'
  'example/ft_realtime_selectiveaverage.m'
  'example/ft_realtime_signalviewer.m'
  'example/ft_realtime_topography.m'
  'online_eeg/ft_realtime_oddball.m'
  'online_eeg/ft_realtime_ouunpod.m'
  'online_eeg/private/compile_midiOut.m'
  'online_eeg/private/controlfunction.m'
  'online_eeg/private/midiOut.m'
  'online_meg/ft_realtime_coillocalizer.m'
  'online_meg/ft_realtime_headlocalizer.m'
  'online_mri/ft_omri_align_init.m'
  'online_mri/ft_omri_align_scan.m'
  'online_mri/ft_omri_info_from_header.m'
  'online_mri/ft_omri_pipeline.m'
  'online_mri/ft_omri_pipeline_nuisance.m'
  'online_mri/ft_omri_quality.m'
  'online_mri/ft_omri_quality_plot.m'
  'online_mri/ft_omri_slice_time_apply.m'
  'online_mri/ft_omri_slice_time_init.m'
  'online_mri/ft_omri_smoothing_kernel.m'
  'online_mri/ft_omri_volume_to_mosaic.m'
  'online_mri/private/coords.m'
  'online_mri/private/encode_nifti1.m'
  'online_mri/private/hom2six.m'
  'online_mri/private/kspace3d.m'
  'online_mri/private/make_A.m'
  'online_mri/private/reslice_vol.m'
  'online_mri/private/rls_init.m'
  'online_mri/private/rls_predict.m'
  'online_mri/private/rls_update.m'
  'online_mri/private/shear_decomp.m'
  'online_mri/private/smooth_vol.m'
  'tutorial/bcifun_latidx.m'
  'tutorial/ft_realtime_asynchronous.m'
  'tutorial/ft_realtime_synchronous.m'
  };

[outlist, depmat] = mydepfun(inlist);

problem = ~cellfun(@isempty, regexp(outlist, 'compat'));
problem = outlist(problem)  % display the output;

if ~isempty(problem)
  error('there are some files that depend on the compat functions');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this was the version of the test script prior to 1 May 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% if isempty(which('ft_realtime_topography'))
%   addpath(fullfile(fileparts(which('ft_defaults')), 'realtime/example'));
% end
%
% cfg = [];
% cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
% cfg.bufferdata = 'first';
% cfg.layout = 'CTF151.lay';
%
% ft_realtime_topography(cfg);


