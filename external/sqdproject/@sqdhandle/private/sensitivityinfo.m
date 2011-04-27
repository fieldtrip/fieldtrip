function t = sensitivityinfo(fid,seekset,channum)
% T = SENSITIVITYINFO(FID,SEEKSET,CHANNUM);
%   FID = valid filepointer to a sqdfile/ sqd-filename
%   SEEKSET = starting point of file read
%   CHANNUM = Channelnumber(s)
% Gets the info regarding senstivity of sensors from the file pointed
% to by fid for given channel_number and returns a structure 

if nargin<1
    error('First argument must be valid file-pointer to a sqd-file');
elseif nargin<2
    seekset = -1; 
end;
filename = '';
if ischar(fid)
    filename = fid;
    fid = fopen(filename,'rb','l');
end;

if fid~=-1
    % Get channelcount 
    fseek( fid, 16, seekset);
    basic_offset = fread(fid,1,'long')+(3*32+8*128+8*128)/8;
    fseek( fid, basic_offset, seekset );
    channelcount	= fread(fid,1,'int');
    if nargin<3
        channum = 0:channelcount-1;
    end;
    if (~isempty(find(channum<0)))|(~isempty(find(channum>channelcount-1)))
        error([ 'Sorry - Invalid channel number, ' num2str(channum) ', given to sensorinfo']);
    end;
    % Get offset of sensitivity values
    fseek( fid, 80, seekset );
    sensitivity_offset = fread(fid,1,'long');
    % Read sensitivity data
    fseek( fid, sensitivity_offset, seekset ); %sizeof(double) = 8 bytes and read 2 doubles
    data = fread(fid,channelcount*2,'double');
    t.Offsets  = [data(channum*2+1)]';
%     t.Gains    = [data(channum*2+2)*10^12]';
    t.Gains    = [data(channum*2+2)*10^15]';
    if ~isempty(filename)
        fclose(fid);
    end;
else
    error('Invalid filename');
end;