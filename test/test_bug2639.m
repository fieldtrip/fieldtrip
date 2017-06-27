function test_bug2639

% TEST ft_checkdata

% MEM 2gb
% WALLTIME 00:10:00

shufflechan = [1 3 2]';

channelcmb = {
  '1' '2'
  '1' '3'
  '2' '3'
  };

%%
freq1o           = [];
freq1o.freq      = 1;
freq1o.label     = {'1', '2', '3'}';
freq1o.dimord    = 'chan_freq';
freq1o.powspctrm = reshape([1 2 3], [3 1]);

freq1r           = freq1o;
freq1r.label     = freq1o.label(shufflechan);
freq1r.powspctrm = freq1o.powspctrm(shufflechan,:);

freq2o = ft_checkdata(freq1o);
freq2r = ft_checkdata(freq1r);
[selo, selr] = match_str(freq2o.label, freq2r.label);
assert(isequal(freq2o.powspctrm(selo,:), freq2r.powspctrm(selr,:)));

%%
freq1o           = [];
freq1o.freq      = 1;
freq1o.cumtapcnt = 1;
freq1o.label     = {'1', '2', '3'}';
freq1o.dimord    = 'rpt_chan_freq';
freq1o.fourierspctrm = reshape([1 2 3], [1 3 1]);

freq1r           = freq1o;
freq1r.label     = freq1o.label(shufflechan);
freq1r.fourierspctrm = freq1o.fourierspctrm(:,shufflechan,:);

freq2o = ft_checkdata(freq1o);
freq2r = ft_checkdata(freq1r);
[selo, selr] = match_str(freq2o.label, freq2r.label);
assert(isequal(freq2o.fourierspctrm(:,selo,:), freq2r.fourierspctrm(:,selr,:)));

%% full, sparse, fourier, sparsewithpow, fullfast

freq2o = ft_checkdata(freq1o, 'cmbrepresentation', 'fourier');
freq2r = ft_checkdata(freq1r, 'cmbrepresentation', 'fourier');
[selo, selr] = match_str(freq2o.label, freq2r.label);
assert(isequal(freq2o.fourierspctrm(:,selo,:), freq2r.fourierspctrm(:,selr,:)));

assert( isequal(freq1r.label,freq2r.label));
assert(~isequal(freq2o.label,freq2r.label));

freq2o = ft_checkdata(freq1o, 'cmbrepresentation', 'full');
freq2r = ft_checkdata(freq1r, 'cmbrepresentation', 'full');
[selo, selr] = match_str(freq2o.label, freq2r.label);
assert(isequal(freq2o.crsspctrm(selo,selo), freq2r.crsspctrm(selr,selr)));

assert( isequal(freq1r.label,freq2r.label));
assert(~isequal(freq2o.label,freq2r.label));

freq2o = ft_checkdata(freq1o, 'cmbrepresentation', 'fullfast');
freq2r = ft_checkdata(freq1r, 'cmbrepresentation', 'fullfast');
[selo, selr] = match_str(freq2o.label, freq2r.label);
assert(isequal(freq2o.crsspctrm(selo,selo), freq2r.crsspctrm(selr,selr)));

assert( isequal(freq1r.label,freq2r.label));
assert(~isequal(freq2o.label,freq2r.label));

freq2o = ft_checkdata(freq1o, 'cmbrepresentation', 'sparse', 'channelcmb', channelcmb);
freq2r = ft_checkdata(freq1r, 'cmbrepresentation', 'sparse', 'channelcmb', channelcmb);
[selo, selr] = match_strcmb(freq2o.labelcmb, freq2r.labelcmb);
assert(isequal(freq2o.crsspctrm(selo,:), freq2r.crsspctrm(selr,:)));

assert( isequal(freq2o.labelcmb,freq2r.labelcmb));

freq2o = ft_checkdata(freq1o, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', channelcmb);
freq2r = ft_checkdata(freq1r, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', channelcmb);
[selo, selr] = match_strcmb(freq2o.labelcmb, freq2r.labelcmb);
assert(isequal(freq2o.crsspctrm(selo,:), freq2r.crsspctrm(selr,:)));
[selo, selr] = match_str(freq2o.label, freq2r.label);
assert(isequal(freq2o.powspctrm(selo,:), freq2r.powspctrm(selr,:)));

assert( isequal(freq1r.label,freq2r.label));
assert(~isequal(freq2o.label,freq2r.label));
assert( isequal(freq2o.labelcmb,freq2r.labelcmb));

%%
freq1o           = [];
freq1o.freq      = 1;
freq1o.cumtapcnt = 1;
freq1o.dimord    = 'chancmb_freq';
freq1o.label     = {'1', '2', '3'}';
freq1o.powspctrm = [1 2 3]';
freq1o.labelcmb  = {
  '1' '2'
  '1' '3'
  '2' '3'
  };
freq1o.crsspctrm = [
  2
  3
  6
  ];

freq1r           = freq1o;
freq1r.label     = freq1o.label(shufflechan);
freq1r.labelcmb  = freq1o.labelcmb(shufflechan,:);
freq1r.powspctrm = freq1o.powspctrm(shufflechan,:);
freq1r.crsspctrm = freq1o.crsspctrm(shufflechan,:);

freq2o = ft_checkdata(freq1o);
freq2r = ft_checkdata(freq1r);
[selo, selr] = match_str(freq2o.label, freq2r.label);
assert(isequal(freq2o.powspctrm(selo,:), freq2r.powspctrm(selr,:)));
[selo, selr] = match_strcmb(freq2o.labelcmb, freq2r.labelcmb);
assert(isequal(freq2o.crsspctrm(selo,:), freq2r.crsspctrm(selr,:)));

freq2o = ft_checkdata(freq1o, 'cmbrepresentation', 'full');
freq2r = ft_checkdata(freq1r, 'cmbrepresentation', 'full');
[selo, selr] = match_str(freq2o.label, freq2r.label);
assert(isequal(freq2o.crsspctrm(selo,selo), freq2r.crsspctrm(selr,selr)));

% assert( isequal(freq1r.label,freq2r.label)); % this is where the channel reordering becomes clear
% assert(~isequal(freq2o.label,freq2r.label)); % this is where the channel reordering becomes clear

freq2o = ft_checkdata(freq1o, 'cmbrepresentation', 'fullfast');
freq2r = ft_checkdata(freq1r, 'cmbrepresentation', 'fullfast');
[selo, selr] = match_str(freq2o.label, freq2r.label);
assert(isequal(freq2o.crsspctrm(selo,selo), freq2r.crsspctrm(selr,selr))); % this is where another problem appears

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sel1, sel2] = match_strcmb(chancmb1, chancmb2)
for i=1:size(chancmb1,1)
  chan1{i} = sprintf('%s_%s', chancmb1{i,:});
end
for i=1:size(chancmb2,1)
  chan2{i} = sprintf('%s_%s', chancmb2{i,:});
end
[sel1, sel2] = match_str(chan1, chan2);

