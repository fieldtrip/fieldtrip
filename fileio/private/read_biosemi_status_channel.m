function [ data ] = read_biosemi_status_channel( statChan )
% This function recieves the BioSemi status channel (the last channel in 
% the data matrix read from the BDF file), and extracts the individual 
% control channels multiplexed within it. See the biosemi website (http://
% www.biosemi.com/faq/trigger_signals.htm) for documentation of this 
% channel structure. 
% 
% The following control channels are read:
% - Trigger Value
% - Extended Trigger Value
% - Start Epoch Flag
% - Speed Value
% - CMS Within Range Flag
% - Battery Low Flag
% - ActiveTwo MK2 Flag
%
% For the trigger channel, the function treats a sequence of idential 
% triggers is read as a single trigger (the first one).
% The function also takes into account the possibility that the status 
% message will change in the middle the sample, resulting in a spurrious 
% intermediate sample where only some of the bits have been updated (for 
% instance, when a trigger value of 10 is sent through a parallel port, the 
% recorded values could be [0 0 0 8 10 10 10], if during the first sample 
% the "2" bit is read before the value was changed but the "8" bit is read 
% after. In such cases, the intermediate value is ignored. Therefore, a
% valid trigger is recorded when: 
% 1) There is a change to a non-zero value followed by zero or more 
% identical values (e.g. the third value in [0 0 10 10 10] or 
% [0 0 10 0 0]). 
% 2) If there are two consecutive non-zero values, the second one is the
% valid trigger (e.g. the third value in [0 8 10 10 10] or [0 8 10 0 0]). 
%
%
% Written by Edden Gerber, lab of Leon Y. Deouell, June 2013


% Define channel bit ranges
bTRIG = 1:8;
bTRIG_EXTENDED = 1:16;
bSTART_EPOCH = 17;
bSPEED = [18:20 22];
bCMS_WITHIN_RANGE = 21;
bLOW_BATT = 23;
bACTIVE_TWO_MK2 = 24;

% Make sure input is column
if isrow(statChan)
    statChan = statChan';
end

% Read trigger channel
v = ExtractBits(statChan,bTRIG);

% Check valid condition flags
val_changed = [v(1)~=0 ; diff(v)~=0];
prev_val_changed = [false ; val_changed(1:end-1)];
val_nonzero = v~=0;
next_val_is_same = [~val_changed(2:end) ; true];
neighbors_are_equal = [ true; v(3:end)-v(1:end-2) == 0 ; true];
next_and_one_before_prev_are_equal = [true ; true ; v(4:end)-v(1:end-3) == 0 ; true];

% Valid condition combinations
valid1 = val_nonzero & val_changed & next_val_is_same; % change to a stable non-zero value
valid2 = val_nonzero & val_changed & neighbors_are_equal; % single-sample change and return to previous value
valid3 = val_nonzero & val_changed & prev_val_changed & next_and_one_before_prev_are_equal; % a single-sample change which is spread across two samples, then returning to the original value

% Either of these is OK
valid = valid1 | valid2 | valid3;

% Delete spurrious triggers
v(~valid) = 0;

% Store in struct
data.Triggers = v;

% Read extended trigger channel (with same steps for removing invalid 
% values)
v = ExtractBits(statChan,bTRIG_EXTENDED);
val_changed = [v(1)~=0 ; diff(v)~=0];
val_nonzero = v~=0;
next_val_is_same = [~val_changed(2:end) ; true];
neighbors_are_equal = [ true; v(3:end)-v(1:end-2) == 0 ; true];
valid1 = val_nonzero & val_changed & next_val_is_same; % change to a stable non-zero value
valid2 = val_nonzero & val_changed & neighbors_are_equal; % single-sample change and return to previous value
valid3 = val_nonzero & val_changed & prev_val_changed & next_and_one_before_prev_are_equal; % a single-sample change which is spread across two samples, then returning to the original value
valid = valid1 | valid2 | valid3;
v(~valid) = 0;
data.ExtTriggers = v;

% Start epoch flag
v = ExtractBits(statChan,bSTART_EPOCH);
change_loc = find(diff([0 ; v]) ~= 0 & diff([v ; 0]) == 0);
t = zeros(size(v));
t(change_loc) = v(change_loc+1);
data.StartEpoch = t;

% Speed
v = ExtractBits(statChan,bSPEED);
change_loc = find(diff([0 ; v]) ~= 0 & diff([v ; 0]) == 0);
t = zeros(size(v));
t(change_loc) = v(change_loc+1);
data.Speed = t;

% CMS is within range flag
v = ExtractBits(statChan,bCMS_WITHIN_RANGE);
change_loc = find(diff([0 ; v]) ~= 0 & diff([v ; 0]) == 0);
t = zeros(size(v));
t(change_loc) = v(change_loc+1);
data.CmsInRange = t;

% Low battery flag
v = ExtractBits(statChan,bLOW_BATT);
change_loc = find(diff([0 ; v]) ~= 0 & diff([v ; 0]) == 0);
t = zeros(size(v));
t(change_loc) = v(change_loc+1);
data.LowBattery = t;

% ActiveTwo MK2 flag
v = ExtractBits(statChan,bACTIVE_TWO_MK2);
change_loc = find(diff([0 ; v]) ~= 0 & diff([v ; 0]) == 0);
t = zeros(size(v));
t(change_loc) = v(change_loc+1);
data.ActiveTwoMk2 = t;

end



function [ out ] = ExtractBits( in, bits)
% This function extracts the specified bits from each number in a vector. 
% 
% in - the input vector
% bits - the bits to be extracted. The least significant bit (LSB) is 1. 
% For example, the decimal 10 is '1010' in binary form, and so 
% ExtractBits(10,[1 2 4]) will extract and concatenate the first, second and fourth bits (counting 
% from the end), resulting in '110' which is 6 in decimal (this is the result of subtracting the 3rd
% bit, correcponding to 4, from 10). 
%
% Written by Edden Gerber, lab of Leon Y. Deouell, June 2013

bits = sort(bits,'ascend');
out = zeros(size(in));
for i=1:length(bits)
    b = floor(mod(in,2^(bits(i)))/2^(bits(i)-1));
    b = b .* 2^(i-1);
    out = out + b;
end
end
