function test_issue1587

% MEM 7gb
% WALLTIME 00:20:00
% DEPENDENCY ft_channelselection

% In https://github.com/fieldtrip/fieldtrip/issues/1587 Alexandra identified an issue
% In https://github.com/fieldtrip/fieldtrip/pull/1599 Jan-Mathijs proposed a solution

% The problem has to do with channel ordering. The code below tests the reading
% functions, which are expected to return the channel order as defined in the file.
% Also the order in hdr.label and grad.label is expected to be the same (although the
% actual channels might be slightly different). The order of the MEG channels in the
% file is non-alphabetic for Neuromag, and alphabetic for CTF.


%%
% ensure that ft_channelselection returns the channel order as in the original data labels

original = {'b', 'a', 'c'}';
selection = 'all';
assert(isequal(ft_channelselection(selection, original), original));
selection = {'a', 'b', 'c'}';
assert(isequal(ft_channelselection(selection, original), original));
selection = {'c', 'b', 'a'}';
assert(isequal(ft_channelselection(selection, original), original));

%%

basedir = '/home/common/matlab/fieldtrip/data/test';

neuromag = {
  './bug1792/20130418_test_cHPI.fif'
  './bug1792/matteo_cHPI.fif'
  './bug1792/sample_chpi.fif'
  './bug1792/sample_chpi_short.fif'
  %  './bug1808/reduced.fif' % this is alphabetic
  %  './bug1914/conversion_testing_MEGEEG_raw.fif'
  %  './bug1914/conversion_testing_MEGandperipheral_raw.fif'
  %  './bug1914/conversion_testing_OnlyMEG_no_triggers_raw.fif'
  %  './bug1914/conversion_testing_OnlyMEG_with_triggers_raw.fif'
  './bug1971/test-ave.fif'
  './bug1971/test_raw.fif'
  './bug2036/sample_audvis_raw.fif'
  './bug2170/an05a4_ss.fif'
  './bug2767/01_ljh_firststd_meg_182.fif'
  './bug731/test_bug731.fif'
  %  './original/meg/neuromag122/jg_single_01raw.fif' % this is alphabetic
  './original/meg/neuromag306/raw.fif'
  './original/meg/neuromag306/run_01_raw.fif'
  './original/rikhenson/Sub15/MEEG/run_01_raw.fif'
  './original/rikhenson/Sub15/MEEG/run_01_sss.fif'
  %  './bug3375/elekta/jn_multimodal_chpi_raw_sss.fif' % this is alphabetic
  };

% neuromag306 datasets are expected to have their channel order non-alphabetical
for i=1:length(neuromag)
  dataset = dccnpath(fullfile(basedir, neuromag{i}));
  hdr = ft_read_header(dataset, 'checkmaxfilter', false);
  grad = ft_read_sens(dataset, 'senstype', 'meg');
  
  meglabel_hdr = ft_channelselection('MEG', hdr.label);
  meglabel_grad = ft_channelselection('MEG', grad.label);
  
  assert(~isequal(hdr.label, sort(hdr.label)));   % not in alphabetical order
  assert(~isequal(grad.label, sort(grad.label))); % not in alphabetical order
  
  [sel1, sel2] = match_str(hdr.label, meglabel_hdr);
  assert(all(diff(sel1)>0)); % monotonously increasing
  assert(all(diff(sel2)>0)); % monotonously increasing
  
  [sel1, sel2] = match_str(meglabel_hdr, meglabel_grad);
  assert(all(diff(sel1)>0)); % monotonously increasing
  assert(all(diff(sel2)>0)); % monotonously increasing
  
end


%%

ctf = {
  './original/meg/ctf151/Subject01.ds'
  './original/meg/ctf275/A0132_Aud-Obj-Recognition_20051115_02.ds'
  './original/meg/ctf275/TacStimRegressConfound.ds'
  %  './original/meg/ctf64/Wat123r1raw.ds' % cannot be read with default implementation
  './bug3280/A1826_comparisonSEF_20161212_02.ds'
  './bug3297/testMEG001_1200hz_20170517_05.ds'
  './bug3375/ctf/muenster_A1331.ds'
  };

% ctf151 and ctf275 datasets are expected to have their channel order alphabetical
for i=1:length(ctf)
  dataset = dccnpath(fullfile(basedir, ctf{i}));
  hdr = ft_read_header(dataset, 'readbids', false);
  grad = ft_read_sens(dataset, 'readbids', false, 'senstype', 'meg');
  
  meglabel_hdr = ft_channelselection('MEG', hdr.label);
  meglabel_grad = ft_channelselection('MEG', grad.label);
  
  assert(isequal(meglabel_hdr, sort(meglabel_hdr)));   % in alphabetical order
  assert(isequal(meglabel_grad, sort(meglabel_grad))); % in alphabetical order
  
  [sel1, sel2] = match_str(hdr.label, meglabel_hdr);
  assert(all(diff(sel1)>0)); % monotonously increasing
  assert(all(diff(sel2)>0)); % monotonously increasing
  
  [sel1, sel2] = match_str(meglabel_hdr, meglabel_grad);
  assert(all(diff(sel1)>0)); % monotonously increasing
  assert(all(diff(sel2)>0)); % monotonously increasing
  
end

%%

dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');

grad = ft_read_sens(dataset, 'senstype', 'meg');

headmodel = [];
headmodel.o = [0 0 40];
headmodel.r = 120;
headmodel.unit = 'mm';
headmodel.type = 'singlesphere';

cfg = [];
cfg.grad = grad;
cfg.headmodel = headmodel;
cfg.sourcemodel.pos = [0 0 70];
cfg.sourcemodel.unit = 'mm';

cfg.channel = 'MEG';
sourcemodel = ft_prepare_leadfield(cfg);
% it should be sorted in alphabetical order (since CTF)
assert(isequal(sourcemodel.label, sort(sourcemodel.label)))

cfg.channel = ft_channelselection('MEG', grad.label);
sourcemodel = ft_prepare_leadfield(cfg);
% it should be sorted in alphabetical order (since CTF)
assert(isequal(sourcemodel.label, sort(sourcemodel.label)))

% The following section does not work as I had expected, since the output is still alphabetical
%
% cfg.channel = flipud(ft_channelselection('MEG', grad.label));
% sourcemodel = ft_prepare_leadfield(cfg);
% % it should NOT be sorted
% assert(~isequal(sourcemodel.label, sort(sourcemodel.label)))

montage = [];
montage.tra = flipud(eye(151));
montage.labelold = ft_channelselection('MEG', grad.label);
montage.labelnew = flipud(ft_channelselection('MEG', grad.label));

cfg.grad = ft_apply_montage(grad, montage);
cfg.channel = 'MEG';
sourcemodel = ft_prepare_leadfield(cfg);
% it should NOT be sorted in alphabetical order, since I flipped it in the grad
assert(~isequal(sourcemodel.label, sort(sourcemodel.label)))
