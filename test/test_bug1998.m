function test_bug1998

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_preprocessing ft_read_data read_neuralynx_ncs

% this bug is detailled on http://bugzilla.fcdonders.nl/show_bug.cgi?id=1998
% and the workaround is explained on http://www.fieldtriptoolbox.org/getting_started/neuralynx?&#discontinuous_recordings

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1998'));

% start with normal preprocessing of a single channel
cfg         = [];
cfg.dataset = 'CSC1.Ncs';
data        = ft_preprocessing(cfg);

% Warning: discontinuous recording, predicted number of timestamps and observed number of timestamps differ by 1693523717.00
%  Please consult the wiki on http://www.fieldtriptoolbox.org/getting_started/neuralynx?&#discontinuous_recordings
% > In fileio/private/read_neuralynx_ncs at 94
%   In ft_read_header at 1196
%   In ft_preprocessing at 394
%
%
% >> disp(data)
%            hdr: [1x1 struct]
%          label: {'CSC1'}
%           time: {[1x12902912 double]}
%          trial: {[1x12902912 double]}
%        fsample: 1893
%     sampleinfo: [1 12902912]
%            cfg: [1x1 struct]

figure
plot(data.time{1})
xlabel('sample number')
ylabel('time (s)')

ts = ft_read_data(cfg.dataset, 'timestamp', 'true'); % raw timestamps
ts = double(ts); % convert from uint64 into double

figure
plot(ts)
xlabel('sample number')
ylabel('timestamp (a.u.)')

% determine the actual timestamps per sample and interpolate the data
dts       = median(diff(ts));
tsinterp  = ts(1):dts:ts(end);
datinterp = interp1(ts, data.trial{1}, tsinterp);

% you can use NaN to replace the data in the gaps
gaps     = find(diff(ts)>2*dts); % skips at least a sample
for igap = 1:length(gaps)
  sel = tsinterp < ts(gaps(igap)+1) & tsinterp > ts(gaps(igap));
  datinterp(:,sel) = NaN;
end

% update the FieldTrip data structure
data.trial{1} = datinterp;
data.time{1}  = (tsinterp - ts(1)) / (data.hdr.Fs .* dts);
