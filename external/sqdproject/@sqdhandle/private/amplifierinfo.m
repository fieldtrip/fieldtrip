function t = amplifierinfo(fid,seekset)
% T = AMPLIFIERINFO(FID,SEEKSET,CHANNUM);
%   FID = valid filepointer to a sqdfile
%   SEEKSET = starting point of file read
%   CHANNUM = Channelnumber
% Gets the info regarding amplier I/O gains from the file pointed
% to by fid for given channel_number and returns a structure

if nargin<1
    error('First argument must be valid file-pointer to a sqd-file');
elseif nargin<2
    seekset = -1;
end;

% Get offset of amplifier information
fseek( fid, 112, seekset); 
amp_offset = fread(fid,1,'long');
% Read amplifier data
t.InputGainBit = -11;
t.InputGainMask= 6144; %(0x1800) in decimal
fseek( fid, amp_offset, seekset );
gain    = fread(fid,1,'int');
% Input gain is stored in Bit-11 to 12
% 0:x1, 1:x2, 2:x5, 3:x10
inputgains = [1,2,5,10];
t.InputGain = inputgains(bitshift(...
    bitand(t.InputGainMask,gain),...
    t.InputGainBit)+1);

% Output gain is stored in Bit-0 to 2
% 0:x1, 1:x2, 2:x5, 3:x10, 4:x20, 5:x50, 6:x100, 7:x200
outputgains = [1,2,5,10,20,50,100,200];
t.OutputGainBit = 0;
t.OutputGainMask= 7;
t.OutputGain = outputgains(bitshift(...
    bitand(t.OutputGainMask,gain),...
    t.OutputGainBit)+1);
