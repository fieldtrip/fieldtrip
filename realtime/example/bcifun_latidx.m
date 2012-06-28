function cmd = bcifun_latidx(cfg,data)
% BCIFUN_LATIDX extracts alpha power, computes the alpha lateralization index 
% and sends a command to an external device
%
% this function computes the lateralization index:
%
% log A / B
%
% for the average alpha power in right channels A and the average alpha
% power in left channels B. In a covert attention experiment we expect
% higher values for attention to the right visual hemifield and lower
% values for attention to the left visual hemifield.
%
% Copyright (C) 2009, Marcel van Gerven

persistent idx;
persistent alis;

if isempty(idx)
    idx = 1; 
    alis = zeros(1,100);
end

% compute channel indices for left and right
rightidx = ismember(data.hdr.label,ft_channelselection({'MR'},data.hdr.label));
leftidx  = ismember(data.hdr.label,ft_channelselection({'ML'},data.hdr.label));

% frequency analysis cfg
cfgf.foi          = 10;
blk               = cfg.blocksize/4;
cfgf.toi          = [blk 2*blk 3*blk];
cfgf.t_ftimwin    = 2*blk;
cfgf.taper        = 'hanning';
 
% compute power and average over time and frequencies
freq  = fastpower(cfgf,data);
right = freq.powspctrm(1,rightidx,:,:);
right = nanmean(right(:));
left  = freq.powspctrm(1,leftidx,:,:);
left = nanmean(left(:));

alis(idx) = log(right/left);
       
if alis(idx) > mean(alis(1:idx)) + 0.5*std(alis(1:idx))
    cmd = 1; % right attention
elseif alis(idx) < mean(alis(1:idx)) - 0.5*std(alis(1:idx))
    cmd = 2; % left attention
else
    cmd = 0; % do nothing
end

if idx==100
    alis = [alis(2:100) 0];
else
    idx = idx+1;
end