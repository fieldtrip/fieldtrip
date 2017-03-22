function test_ft_datatype

% MEM 8gb
% WALLTIME 02:00:00

% TEST ft_datatype ft_datatype_comp ft_datatype_mvar ft_datatype_source ft_datatype_dip ft_datatype_parcellation ft_datatype_spike ft_datatype_freq ft_datatype_raw ft_datatype_timelock ft_datatype_headmodel ft_datatype_segmentation ft_datatype_volume ft_datatype ft_datatype_sens

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this style is also used in test_ft_analysisprotocol and test_ft_datatype_source

dirlist = {
  dccnpath('/home/common/matlab/fieldtrip/data/test/latest')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20131231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20130630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20121231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20120630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20111231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20110630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20101231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20100630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20091231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20090630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20081231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20080630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20071231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20070630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20061231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20060630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20051231')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20050630')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20040623')
  dccnpath('/home/common/matlab/fieldtrip/data/test/20031128')
  };

for j=1:length(dirlist)
  filelist = hcp_filelist(dirlist{j});
  
  [dummy, dummy, x] = cellfun(@fileparts, filelist, 'uniformoutput', false);
  sel = strcmp(x, '.mat');
  filelist = filelist(sel);
  clear p f x
  
  nbytes = zeros(numel(filelist),1);
  for i=1:length(filelist)
    % sort for the file size
    tmp = dir(filelist{i});
    nbytes(i) = tmp.bytes;
  end
  [nbytes,ix] = sort(nbytes);
  filelist = filelist(ix);
  
   for i=1:length(filelist)
    
    try
      fprintf('processing data structure %d from %d\n', i, length(filelist));
      fprintf('file size %d\n', nbytes(i));
      var = loadvar(filelist{i});
    catch
      % some of the mat files are corrupt, this should not spoil the test
      % disp(var);
      disp(lasterr);
      continue
    end
    
    type = 'unknown';
    
    if ~isempty(regexp(filelist{i}, '/raw/'))
      type = 'raw';
    elseif ~isempty(regexp(filelist{i}, '/comp/'))
      type = 'comp';
    elseif ~isempty(regexp(filelist{i}, '/timelock/'))
      type = 'timelock';
    elseif ~isempty(regexp(filelist{i}, '/freq/'))
      type = 'freq';
    elseif ~isempty(regexp(filelist{i}, '/source/'))
      type = 'source';
    elseif ~isempty(regexp(filelist{i}, '/volume/')) || ~isempty(regexp(filelist{i}, '/mri/'))
      type = 'volume';
    end
    
    switch type
      case 'raw'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      case 'comp'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
        type = 'raw'; % comp data is a special type of raw data
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      case 'timelock'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      case 'freq'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      case 'source'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      case 'volume'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      otherwise
        warning('not testing %s', filelist{i});
        % do nothing
    end % switch
    clear var type
    
  end % for filelist
end % for dirlist


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see how it deals with raw/timelock/freq structures with component topographies
% see SVN revision 9336 and http://bugzilla.fcdonders.nl/show_bug.cgi?id=2518

raw = [];
raw.label = {'1', '2', '3'};
for i=1:5
  raw.time{i} = 1:10;
  raw.trial{i} = randn(3, 10);
end

timelock = [];
timelock.label = {'1', '2', '3'};
timelock.time = 1:10;
timelock.avg = randn(3, 10);
timelock.dimord = 'chan_time';

freq = [];
freq.label = {'1', '2', '3'};
freq.freq = 1:5;
freq.time = 1:10;
freq.powspctrm = randn(3, 5, 10);
freq.dimord = 'chan_freq_time';

assert(ft_datatype(raw,       'raw'), 'incorrect datatype');
assert(ft_datatype(timelock,  'timelock'), 'incorrect datatype');
assert(ft_datatype(freq,      'freq'), 'incorrect datatype');

% create some structures that have features of two different types

rawc = raw;
rawc.topo = randn(4,3);
rawc.unmixing = randn(3,4);
rawc.topolabel = {'a', 'b', 'c', 'd'};

timelockc = timelock;
timelockc.topo = randn(4,3);
timelockc.unmixing = randn(3,4);
timelockc.topolabel = {'a', 'b', 'c', 'd'};

freqc = freq;
freqc.topo = randn(4,3);
freqc.unmixing = randn(3,4);
freqc.topolabel = {'a', 'b', 'c', 'd'};

assert(ft_datatype(rawc,       'raw'),      'incorrect datatype');
assert(ft_datatype(timelockc,  'timelock'), 'incorrect datatype');
assert(ft_datatype(freqc,      'freq'),     'incorrect datatype');

assert(ft_datatype(rawc,       'comp'), 'incorrect datatype');
assert(ft_datatype(timelockc,  'comp'), 'incorrect datatype');
assert(ft_datatype(freqc,      'comp'), 'incorrect datatype');

% the default is that for raw+comp ft_datatype returns comp
% the default for timelock+comp and freq+comp is to return them as timelock or freq for backward compatibility
assert(strcmp(ft_datatype(rawc       ), 'raw+comp'),      'raw+comp datatype was expected');
assert(strcmp(ft_datatype(timelockc  ), 'timelock+comp'), 'timelock+comp datatype was expected');
assert(strcmp(ft_datatype(freqc      ), 'freq+comp'),     'freq+comp datatype was expected');

