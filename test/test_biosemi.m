% function test_biosemi

% DEPENDENCY read_biosemi_bdf write_biosemi_bdf
% WALLTIME 00:10:00
% MEM 2gb
% DATA private

%%

prefix = dccnpath('/project/3031000.02/test/original/eeg/bdf');

filelist = {
  '050327BH_overCZnoAlpha.bdf'
  'mbrain-train-smarting-mobile-eeg.bdf'
  'Newtest17-2048.bdf'
  'Newtest17-256.bdf'
  %'test_generator_2.bdf' % this has channels with different sampling rates
  };

formatlist = {
  'biosemi_v1'
  'biosemi_v2'
  'biosemi_v3'
  ''            % the default, whatever that is
  'biosig'      % this should also work
  };

clear hdr dat event

for j=1:numel(formatlist)
  format = formatlist{j};
  for i=1:numel(filelist)
    filename = fullfile(prefix, filelist{i});
    fprintf('reading file %s as %s\n', filename, format);
    hdr{j,i}   = ft_read_header(filename, 'headerformat', format);
    dat{j,i}   = ft_read_data(filename, 'dataformat', format);
    event{j,i} = ft_read_event(filename, 'eventformat', format);
  end
end

%%

for i=1:numel(filelist)
  for j=1:3 % only compare the biosemi_vN
    disp([i, j])
    assert(hdr{j,i}.Fs==hdr{1,i}.Fs);
    assert(hdr{j,i}.nChans==hdr{1,i}.nChans);
    assert(hdr{j,i}.nSamples*hdr{j,i}.nTrials==hdr{1,i}.nSamples*hdr{1,i}.nTrials);
  end
end

%%

tempfile = [tempname '.bdf'];

for i=1:numel(filelist)
  for j=1:numel(formatlist)
    disp([i, j])

    ft_write_data(tempfile, dat{j,i}, 'header', hdr{j,i}, 'dataformat', 'biosemi_bdf');
    tmphdr   = ft_read_header(tempfile, 'headerformat', 'biosemi_v3');
    tmpdat   = ft_read_data(tempfile, 'dataformat', 'biosemi_v3');
    tmpevent = ft_read_event(tempfile, 'eventformat', 'biosemi_v3');

    assert(tmphdr.Fs==hdr{j,i}.Fs);
    assert(tmphdr.nChans==hdr{j,i}.nChans);
    assert(tmphdr.nSamples*tmphdr.nTrials==hdr{j,i}.nSamples*hdr{j,i}.nTrials);

  end
end

% clean up the temporary file
delete(tempfile);
