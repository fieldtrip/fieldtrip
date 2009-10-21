function [hdr, dat] = read_mpi_ds(dirname)

% READ_MPI_DS reads all DAP files from a directory containing files or
% alternatively a single DAP file and returns it in a simplified FieldTrip
% format. The analog channels and spike channels are both returned in a
% continuous format.
%
% Use as
%   [hdr, dat] = read_mpi_ds(dirname)
% or
%   [hdr, dat] = read_mpi_ds(filename)
%
% See also READ_MPI_DAP

% Copyright (C) 2005-2007, Robert Oostenveld
%
% $Log: read_mpi_ds.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.4  2007/03/19 16:53:49  roboos
% allow for multiple spikes at the same sample (i.e. increment by 1 for each spike)
%

hdr = [];
dat = [];

if isdir(dirname)
  ls = dir(dirname);
  dapfile = {};
  for i=1:length(ls)
    if ~isempty(regexp(ls(i).name, '.dap$'))
      dapfile{end+1} = fullfile(dirname, ls(i).name);
    end
  end
  dapfile = sort(dapfile);
elseif iscell(dirname)
  dapfile = dirname;
else
  dapfile = {dirname};
end

for i=1:length(dapfile)
  dap        = read_mpi_dap(dapfile{i});
  siz        = size(dap.analoghdr);
  nanalog    = dap.filehdr.nanalog;
  nspike     = dap.filehdr.nspike;
  nexp(i)    = siz(1);
  nsweep(i)  = siz(2);
  nsample(i) = length(dap.analog{1});
  analog     = zeros(nanalog, nsample(i)*nsweep(i)*nexp(i));
  spike      = zeros(nspike , nsample(i)*nsweep(i)*nexp(i));
  % only collect the data when required
  if nargout>1
    % collect all analog data of all sweeps in a single matrix
    for e=1:nexp(i)
      for s=1:nsweep(i)
        for c=1:nanalog
          begsample = (e*s-1) * nsample(i) + 1;
          endsample = (e*s  ) * nsample(i);
          analog(c,begsample:endsample) = dap.analog{e,s,c}(:)';
        end
      end
    end
    % collect all spike data of all sweeps in a single matrix
    for e=1:nexp(i)
      for s=1:nsweep(i)
        for c=1:nspike
          spike(c,dap.spike{e,s,c} + (e*s-1) * nsample(i)) = spike(c,dap.spike{e,s,c} + (e*s-1) * nsample(i)) + 1;
        end
      end
    end
    dat = cat(1, analog, spike);
  end % if nargout
end

hdr.Fs          = 1000;              % sampling frequency
hdr.nChans      = nanalog+nspike;    % number of channels
hdr.nSamples    = nsample(1);        % number of samples per trial
hdr.nSamplesPre = 0;                 % number of pre-trigger samples in each trial
hdr.nTrials     = sum(nsweep.*nexp); % number of trials
hdr.label       = {};                % cell-array with labels of each channel
for i=1:nanalog
  hdr.label{i} = sprintf('analog%02d', i);
end
for i=1:nspike
  hdr.label{nanalog+i} = sprintf('spike%02d', i);
end

