function cmd = bcifun_latidx(cfg,data)

% BCIFUN_LATIDX extracts alpha power, computes the alpha lateralization index 
% and sends a command to an external device
%
% This function computes the lateralization index log(A/B) for the
% average alpha power in right channels A and the average alpha power
% in left channels B. In a covert attention experiment we expect
% higher values for attention to the right visual hemifield and lower
% values for attention to the left visual hemifield.

% Copyright (C) 2009, Marcel van Gerven
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

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
