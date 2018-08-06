function data_tstamps = comp_tstamps(inp, sfreq)
% COMP_TSTAMPS - extract timestamps from a trigger channel
%   INP - vector of samples for the trigger channel
%   SFREQ - sampling frequency
%   Return the vector of the same length as INP, containing timestamps for
%   each entry of INP. For detecting timestamps use parameters defined
%   below (should match the parameters used for generating the timing
%   sequence).
%
%   TODO: this function does not handle the boundary case for the first train
%   of pulses correctly. This is because there is no trigger before the train
%   and there will be no dtrigs value before the first trigger of the train.
%   Thus the first pulse train will always be ignored. It would be neat to fix
%   this.

%--------------------------------------------------------------------------
%   Copyright (C) 2015 BioMag Laboratory, Helsinki University Central Hospital
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, version 3.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------

THRESH = 3;
BASELINE = 5;   % seconds
TRAIN_STEP = 0.015; % seconds
NBITS = 43;          % including the parity bit

% find the upgoing flank of all triggers (threshold crossings)
trigs = find([0 diff(inp>THRESH)>0]);

% determine the duration of each trigger
d = trigs(2:end) - trigs(1:end-1);

samps = [];
tss = [];

% iterate over all timestamp candidates
for i = find(d > BASELINE * sfreq)
    ts = read_timestamp(d, i, TRAIN_STEP*sfreq, NBITS);
    if ts ~= -1
        samps(end+1) = trigs(i+1);
        tss(end+1) = ts;
    end
end
    
% fit timestamps to samples with linear regression
p = polyfit(samps, tss, 1);
data_tstamps = [1:length(inp)] * p(1) + p(2);

% DEBUG
%figure;
%title('Fit accuracy');
%plot(samps/sfreq, p(1)*samps+p(2) - tss);
%xlabel('time, seconds');
%ylabel('linear fit error, msec');
% ~DEBUG


function ts = read_timestamp(dtrigs, cur, step, nbits)
% READ_TIMESTAMP - read and decode one timestamp

ts = 0;
parity = false;

for i = 1 : nbits
    % end of input reached before NBITS bits read
    if cur+i > length(dtrigs)
        warning('end of input reached before NBITS bits read');
        ts = -1;
        return;
    end
    
    % invalid interval between two triggers
    if (dtrigs(cur+i) < step*1.5) | (dtrigs(cur+i) > step*4.5)
        warning('invalid interval between two triggers');
        ts = -1;
        return;
    end
    
    if dtrigs(cur+i) > step*3
        parity = ~parity;
        if i < nbits    % don't read the parity bit into the timestamp
            ts = ts + 2^(i-1);
        end
    end
end

if parity
    warning('parity check failed');
    ts = -1;
end
