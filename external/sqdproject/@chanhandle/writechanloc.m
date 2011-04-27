function writechanloc(t,locfile,locform)
% WRITECHANLOC(CHANHANDLE_OBJ,FILENAME); Writes the location of the channels in the
% CHANHANDLE object to FILENAME.
%
% WRITECHANLOC(CHANHANDLE_OBJ,FILENAME,FORMAT);
% where FORMAT can be:
%   'xyz' : (Default) Matlab/EEGLAB Cartesian coordinates (NOT EGI Cartesian).
%           x is toward nose; y is toward left ear; z is toward vertex
%
% This is useful for analysing the data in EEGLAB toolbox
% Please see  http://sccn.ucsd.edu/eeglab for more details of the toolbox.

if nargin<1
    error('Atleast one input (chanhandle object) required');
end;

if ~isa(t,'chanhandle')
    error('First input must belong to CHANHANDLE class');
end;

switch nargin   % init input arguments
case 1
    locfile = 'testlocation.xyz';
    locform = 'xyz';
case 2
    locform = 'xyz';
end;

switch locform
case 'xyz'      % if output format = 'xyz' ie Matlab/EEGLAB Cartesian coordinates
    locstr = '';
    try
        % read requird info (x,y,z locations from object)
        % location might not be specified for some channels
        % hence, the try-catch
        chanlocx    = get(t,'SensorInfo.x');
        chanlocy    = get(t,'SensorInfo.y');
        chanlocz    = get(t,'SensorInfo.z');
    catch
        disp('Sorry - Write Failed');
        disp('Not all of the given Channels have positions specified');
        error('Please check Channel Numbers to verify');
    end;
    channum     = get(t,'ChannelNumber');   % get hardware channel number
    chanlabel   = get(t,'SensorInfo.type'); % get sensor type
    
    for i = 1:get(t,'Numchannels')  % Write to output string for each channel
        locstr = [locstr,sprintf('%d\t\t%1.3f\t\t%1.3f\t\t%1.3f\t\t%s\n',...
                channum(i)+1,chanlocx(i),chanlocy(i),chanlocz(i),chanlabel{i})];
    end;
    fid = fopen(locfile,'w','b');   % Open as bigendian
    if fid~=-1
        fwrite(fid,locstr);         % Write string to file
        fclose(fid);
    else
        error('Write Unsuccessful - File not found');
    end;
otherwise
    error(['Sorry - ''' locform ''' not supported yet. See help WRITECHANLOC']);
end;