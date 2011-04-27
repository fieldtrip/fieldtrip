function t = sensorinfo(fid,seekset,channum)
% T = SENSORINFO(FID,SEEKSET,CHANNUM);
%   FID = valid filepointer to a sqdfile
%   SEEKSET = starting point of file read
%   CHANNUM = CHANNEL_NUMBER
% Gets the info regarding sensor location from the file pointed
% to by fid for the given channel number and returns a structure

if nargin<1
    error('First argument must be valid file-pointer to a sqd-file');
elseif nargin<2
    channum = 0;
    seekset = -1;
elseif nargin<3
    channum = 0;
end;


%********** Channel Information **********
% Get channelcount 
fseek( fid, 16, seekset);
basic_offset = fread(fid,1,'long')+(3*32+8*128+8*128)/8;
fseek( fid, basic_offset, seekset );
channelcount	= fread(fid,1,'int');
if (channum<0)|(channum>channelcount-1)
    error([ 'Sorry - Invalid channel number, ' num2str(channum) ', given to sensorinfo']);
end;
% Get offset of channel information
fseek( fid, 64, seekset );
offset   = fread(fid,1,'long');
size     = fread(fid,1,'long');

% Read and Output channel info
fseek(fid,offset+size*channum,seekset);
type = fread(fid,1,'int');
switch type
case 1  % magnetometer
    t.type = 'magnetometer';
    t.x = fread(fid,1,'double');
    t.x = fread(fid,1,'double');
    t.z = fread(fid,1,'double');
    t.zdir = fread(fid,1,'double');
    t.xdir = fread(fid,1,'double');
    t.coilsize = fread(fid,1,'double');
case 2  % axialgradiometer
    t.type = 'axialgradiometer';
    t.x = fread(fid,1,'double');
    t.y = fread(fid,1,'double');
    t.z = fread(fid,1,'double');
    t.zdir = fread(fid,1,'double');
    t.xdir = fread(fid,1,'double');
    t.baseline = fread(fid,1,'double');
    t.coilsize = fread(fid,1,'double');
case 3  % plannergradiometer
    t.type = 'plannergradiometer';
    t.x = fread(fid,1,'double');
    t.y = fread(fid,1,'double');
    t.z = fread(fid,1,'double');
    t.zdir = fread(fid,1,'double');
    t.xdir = fread(fid,1,'double');
    t.zdir2 = fread(fid,1,'double');
    t.xdir2 = fread(fid,1,'double');
    t.baseline = fread(fid,1,'double');
    t.coilsize = fread(fid,1,'double');
case 4  % magnetometer(reference)
    t.type = 'magnetometer_reference'
    t.x = fread(fid,1,'double');
    t.y = fread(fid,1,'double');
    t.z = fread(fid,1,'double');
    t.zdir = fread(fid,1,'double');
    t.xdir = fread(fid,1,'double');
    t.coilsize = fread(fid,1,'double');
case 5  % axialgradiometer(reference)
    t.type = 'axialgradiometer_reference';
    t.x = fread(fid,1,'double');
    t.y = fread(fid,1,'double');
    t.z = fread(fid,1,'double');
    t.zdir = fread(fid,1,'double');
    t.xdir = fread(fid,1,'double');
    t.baseline = fread(fid,1,'double');
    t.coilsize = fread(fid,1,'double');
case 6  % planergradiometer(reference)
    t.type = 'planergradiometer_reference';
    t.x = fread(fid,1,'double');
    t.y = fread(fid,1,'double');
    t.z = fread(fid,1,'double');
    t.zdir = fread(fid,1,'double');
    t.xdir = fread(fid,1,'double');
    t.zdir2 = fread(fid,1,'double');
    t.xdir2 = fread(fid,1,'double');
    t.baseline = fread(fid,1,'double');
    t.coilsize = fread(fid,1,'double');
case -1 % trigger
    t.type = 'trigger';
    t.identnum = fread(fid,1,'int');
    t.channame = fread(fid,64/8,'uchar');
    t.channame = setstr(t.channame');
case -2 % EEG
    t.type = 'EEG';
    t.identnum = fread(fid,1,'int');
    t.channame = fread(fid,64/8,'uchar');
    t.channame = setstr(t.channame');
case -3     % ECG
    t.type = 'ECG';
    t.identnum = fread(fid,1,'int');
    t.channame = fread(fid,64/8,'uchar');
    t.channame = setstr(t.channame');
case -4     % etc.
    t.type = 'etc';
    t.identnum = fread(fid,1,'int');
    t.channame = fread(fid,64/8,'uchar');
    t.channame = setstr(t.channame');
otherwise   % NULL
    t.type = 'NULL';
    t.identnum = '';
    t.channame = 'NULL';
end;    % End switch on type
