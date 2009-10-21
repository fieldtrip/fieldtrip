function [stimulus, response, segment, timezero] = read_brainvision_vmrk(filename);

% READ_BRAINVISION_VMRK reads the markers and latencies
% it returns the stimulus/response code and latency in ms.
%
% Use as
%   [stim, resp, segment, timezero] = read_brainvision_vmrk(filename)
% 
% This function needs to read the header from a separate file and
% assumes that it is located at the same location.
%
% See also READ_BRAINVISION_VHDR, READ_BRAINVISION_EEG

% original M. Schulte 31.07.2003
% modifications R. Oostenveld 14.08.2003
%
% $Log: read_brainvision_vmrk.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.3  2008/07/24 12:05:12  roboos
% changed end ot line to unix style
%
% Revision 1.2  2004/03/30 07:23:28  roberto
% added to CVS repository, copyrights added and implemented multiple file
% formats
%

stimulus=[];
response=[];
segment=[];
timezero=[];

% read the header belonging to this marker file
hdr=read_brainvision_vhdr([filename(1:(end-4)) 'vhdr']);

fid=fopen(filename,'rt');
if fid==-1,
    error('cannot open marker file')
end

line=1;
while line~=-1,
  line=fgetl(fid);
  % pause
  if ~isempty(line),
    if ~isempty(findstr(line,'Mk')),
      if ~isempty(findstr(line,'Stimulus'))
        [token,rem] = strtok(line,',');
        type=sscanf(rem,',S %i');
        [token,rem] = strtok(rem,',');
        time=(sscanf(rem,', %i')-1)/hdr.Fs*1000; 
        stimulus=[stimulus; type time(1)];

      elseif ~isempty(findstr(line,'Response'))
        [token,rem] = strtok(line,',');
        type=sscanf(rem,',R %i');
        [token,rem] = strtok(rem,',');
        time=(sscanf(rem,', %i')-1)/hdr.Fs*1000;                
        response=[response; type, time(1)];

      elseif ~isempty(findstr(line,'New Segment'))
        [token,rem] = strtok(line,',');
        time=(sscanf(rem,',,%i')-1)/hdr.Fs*1000;                
        segment=[segment; time(1)];

      elseif ~isempty(findstr(line,'Time 0'))
        [token,rem] = strtok(line,',');
        time=(sscanf(rem,',,%i')-1)/hdr.Fs*1000;                
        timezero=[timezero; time(1)];

      end
    end
  else 
    line=1;
  end
end

fclose(fid);    

