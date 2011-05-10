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
fseek( fid, amp_offset, seekset );
gain    = fread(fid,1,'int');

% 2 amplifier stage / 12bit ADC systems
if gain < hex2dec('FFFF')
    t.InputGainBit = -11;
    t.InputGainMask= 6144; %(0x1800) in decimal
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
    
% 3 amplifier stage / 16bit ADC systems
else
    % gain
    % 0:x1, 1:x2, 2:x5, 3:x10, 4:x20, 5:x50, 6:x100, 7:x200
    gains = [1,2,5,10,20,50,100,200];
    % gain1 is stored in Bit-12 to 15
    Gain1Bit = -12;
    Gain1Mask= hex2dec('00007000');
    Gain1 = gains(bitshift(...
        bitand(Gain1Mask,gain),...
        Gain1Bit)+1);
   
    t.InputGain = Gain1;
     
    % gain2 is stored in Bit 28 to 31
    Gain2Bit = -28;
    Gain2Mask= hex2dec('70000000');
    Gain2 = gains(bitshift(...
        bitand(Gain2Mask,gain),...
        Gain2Bit)+1);
   
  
    % gain3 is stored in Bit 24 to 27
    Gain3Bit = -24;
    Gain3Mask= hex2dec('07000000');
    Gain3 = gains(bitshift(...
        bitand(Gain3Mask,gain),...
        Gain3Bit)+1);
    
    % take care of the ratio 12bit ADC to the new 16bit ADC
    % 12bit is hard coded in line 130 of getdata.m
    t.OutputGain = Gain2*Gain3 *(2^16)/(2^12);
end