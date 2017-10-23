function test_yokogawa

% MEM 1500mb
% WALLTIME 00:10:00

% TEST hasyokogawa read_yokogawa_data read_yokogawa_event read_yokogawa_header yokogawa2grad yokogawa2vol

% this script tests some files from the three different types of yokogawa MEG systems
% it tests the general reading and whether the system type and channel selection all work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/yokogawa160/Continuous1.con');
hdr = ft_read_header(filename);

if length(ft_channelselection('MEG', hdr.label))~=160
  error('did not select all MEG channels');
end

if ~ft_senstype(hdr, 'yokogawa')
  error('hdr ~= yokogawa');
elseif ~ft_senstype(hdr.grad, 'yokogawa')
  error('grad ~= yokogawa');
elseif ~ft_senstype(hdr.label, 'yokogawa')
  error('label ~= yokogawa');
end

if ~ft_senstype(hdr, 'yokogawa160')
  error('hdr ~= yokogawa160');
elseif ~ft_senstype(hdr.grad, 'yokogawa160')
  error('grad ~= yokogawa160');
elseif ~ft_senstype(hdr.label, 'yokogawa160')
  error('label ~= yokogawa160');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/yokogawa440/S1_MEG_Epoch.raw');
hdr = ft_read_header(filename);

if length(ft_channelselection('MEG', hdr.label))<400
  % it should have at least 400 MEG channels
  error('did not select all MEG channels');
end

if ~ft_senstype(hdr, 'yokogawa')
  error('hdr ~= yokogawa');
elseif ~ft_senstype(hdr.grad, 'yokogawa')
  error('grad ~= yokogawa');
elseif ~ft_senstype(hdr.label, 'yokogawa')
  error('label ~= yokogawa');
end

if ~ft_senstype(hdr, 'yokogawa440')
  error('hdr ~= yokogawa440');
elseif ~ft_senstype(hdr.grad, 'yokogawa440')
  error('grad ~= yokogawa440');
elseif ~ft_senstype(hdr.label, 'yokogawa440')
  error('label ~= yokogawa440');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/yokogawa64/2011_01_28_0354_ME053_AEF.con');
hdr = ft_read_header(filename);

if length(ft_channelselection('MEG', hdr.label))~=64
  error('did not select all MEG channels');
end

if ~ft_senstype(hdr, 'yokogawa')
  error('hdr ~= yokogawa');
elseif ~ft_senstype(hdr.grad, 'yokogawa')
  error('grad ~= yokogawa');
elseif ~ft_senstype(hdr.label, 'yokogawa')
  error('label ~= yokogawa');
end

if ~ft_senstype(hdr, 'yokogawa64')
  error('hdr ~= yokogawa64');
elseif ~ft_senstype(hdr.grad, 'yokogawa64')
  error('grad ~= yokogawa64');
elseif ~ft_senstype(hdr.label, 'yokogawa64')
  error('label ~= yokogawa64');
end

