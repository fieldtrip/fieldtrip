function trl = ft_trialfun_neuralynx(cfg)

% FT_TRIALFUN_NEURALYNX outputs trials for Neuralynx datasets.
%
% Use as
%   cfg          = [];
%   cfg.dataset  = 'directory_with_ncs_files';
%   cfg.trialfun = 'ft_trialfun_neuralynx';
%   cfg = ft_definetrial(cfg);
%   data = ft_preprocessing(cfg);
%
% A Neuralynx dataset consists of a directory containing as many .ncs files
% as electrodes. The option cfg.dataset should refer to this directory.
%
% Neuralynx data is saved as a series of 512-sample blocks (or "records").
% When user presses pause or stop buttons, the recording stops at sample i
% from the current record, and the remaining samples (i+1 to 512) are a
% copy of samples i+1 to 512 of the previous record. These copied samples
% are marked as invalid in ncs.NumValidSamples, and have to be discarded as
% they are not genuine samples and as they introduce redundancies.
% 
% Contact ludovic.bellier@berkeley.edu for any questions.


gapThreshold = 3; % in number of records (e.g., at fs=8kHz, 3*512 samples represents 192 ms)
recordsize = 512;

dirname = cfg.dataset;
if ~any(strcmp(dirname(end), {'/', '\'}))
    dirname = [dirname filesep];
end

ls = dir([dirname '*ncs']);
fname = cellfun(@(x) fullfile(dirname, x), {ls(:).name}, 'un', 0);

ncs = read_neuralynx_ncs(fname{1});
idxPauses = find(ncs.NumValidSamp < recordsize);
nValidSamples = ncs.NumValidSamp(idxPauses);
nPauses = length(idxPauses);
isPause = false(nPauses, 1);
for i = 1:nPauses
    if idxPauses(i) == ncs.NRecords
        isPause(i) = true;
    else
        timeGap = double(diff(ncs.TimeStamp(idxPauses(i):idxPauses(i)+1)));
        if timeGap > gapThreshold * mode(diff(ncs.TimeStamp))
            isPause(i) = true;
        end
    end
end
idxPauses = idxPauses(isPause);
nValidSamples = nValidSamples(isPause);
nTrials = length(idxPauses);
trl = zeros(nTrials, 3);
trl(:, 2) = ((idxPauses-1)*recordsize + nValidSamples)';
trl(1, 1) = 1;
trl(2:end, 1) = (idxPauses(1:end-1)*recordsize+1)';