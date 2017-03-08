function test_bug1856

% MEM 1500mb
% WALLTIME 00:20:00

% TEST ft_read_header ft_read_sens ft_chanunits

% make sure that grad.chanunit and hdr.chanunit are specified for the most
% common MEG systems, biosemi, brainvision and egi

% FIXME egi example data is missing in this list

filename = {
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/bdf/Newtest17-256.bdf')
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/brainvision/Mischa.vhdr')
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/bti148/c,rfhp0.1Hz.m4d')
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/bti248/e,rfDC')
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/bti248grad/e,rfhp1.0Hz,COH')
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf151/Subject01.ds')
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf275/A0132_Aud-Obj-Recognition_20051115_02.ds')
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag122/jg_single_01raw.fif')
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/raw.fif')
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/yokogawa160/Continuous1.con')    % this one has 28 channels of an unknown type
  dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/bdf/050327BH_overCZnoAlpha.bdf') % this has 10-20 EEG channel labels, see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1856#c8
  };

for i=1:length(filename)
  hdr = ft_read_header(filename{i});

  assert(isfield(hdr, 'chantype'));
  assert(isfield(hdr, 'chanunit'));
  unknown = strcmp(hdr.chanunit, 'unknown');
  if mean(unknown)>0.30
    % allow for 20% unknown
    error('there are too many channels with unknown units in the header structure');
  end

  if ~isempty(regexp(filename{i}, 'meg', 'once'))
    grad = ft_read_sens(filename{i}, 'senstype', 'meg'); % ensure that the grad is returned for the 306 neuromag file
    assert(isfield(grad, 'chantype'));
    assert(isfield(grad, 'chanunit'));
    unknown = strcmp(grad.chanunit, 'unknown');
    if any(unknown)
      % do not allow any unknown
      error('there are too many channels with unknown units in the grad structure');
    end
  end % if meg

end

