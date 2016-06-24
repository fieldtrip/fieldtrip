function [out] = read_ced_son(datafile,varargin)

% READ_CED_SON
%
% [OUT] = read_ced_son(DATAFILE,VARARGIN);
%
% Reads a analog and event data from a CED SON file 
% (SON files are created by Spike2 software). Currently, only 
% analog channels and event data can be read. 
%
% Optional parameter      Default
%   'readevents'          'no'
%   'readdata'            'no'
%   'readtimestamps'      'no'
%   'begsample'           -1
%   'endsample'           -1
%   'channels'             []      
% 
% Please note that CED DAQ systems do a sequential ADC, thus 
% channels do not share the same time axis: The timestamps of the
% analog channels differ on a subsample level. Use the 'readtimestamps'
% input parameter to get a matrix with time axes corresponding
% to the data channels.
%
% Use begsample and endsample parameters to specify the boundaries
% of the requested data chunk. Setting these parameters to -1 will
% return data from the start or until the end of the datafile,
% respectively.
%
% Specifying [1,2] for 'channels' will load the 1st and the 2nd 
% analog channel, __regardless of the actual channel number__
% If, for example channel 1,2,3 are event channels, 4 as an analog
% channel, 5 is an event channel, and 6 is and analog channel, 
% specifying [1 2] for 'channels' will load analog channel 4 and 6.
% Specifying [] for channels will return all analog channels.
%
% Setting 'readtimestamps' to 'yes' will return a time vector for
% each analog channel.
%
% Depending on the input parameters, the function will return a structure 
% with fields: 
%    'header'       Header information of the SON file
%    'event'        All data from event channels are pooled 
%                   and stored in this structure.
%    'data'         Cell-array with analog data
%    'time'         Cell-array with time vectors corresponding to 'data'
%
% Uses Neuroshare libraries to read Spike2 SON data 
% (see: http://neuroshare.sourceforge.net)

% Gijs van Elswijk - 2005 (v0.1)

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

MODE = 'continuous';  % assume continuous now

% process in put parameters
if ~isempty(varargin) pars=struct(varargin{:}); else pars=[]; end;
if ~isfield(pars,'readdata'),       pars = setfield(pars,'readdata','no');       end;
if ~isfield(pars,'readevents'),     pars = setfield(pars,'readevents','no');     end;
if ~isfield(pars,'readtimestamps'), pars = setfield(pars,'readtimestamps','no');  end;
if ~isfield(pars,'begsample'),      pars = setfield(pars,'begsample',-1);        end;
if ~isfield(pars,'endsample'),      pars = setfield(pars,'endsample',-1);        end;
if ~isfield(pars,'channels'),       pars = setfield(pars,'channels',[]);         end;

% set all fields string values to lowercase
fields = fieldnames(pars);
for idx=1:length(fields)
    if ischar(getfield(pars,fields{idx})),
        pars=setfield(pars,fields{idx},lower(getfield(pars,fields{idx})));
    end;
end;

% First, check if NeuroShare DLL can be loaded
if ft_filetype(datafile, 'ced_son')
    % TODO other DLLs for other binary formats could be supported here as well
    ns_RESULT = ns_SetLibrary(which('nsCedSon.dll'));
end
[st,libinfo] = ns_GetLibraryInfo;
if st,
    error(['Could not get NeuroShare library info, please use the NS_SETLIBRARY function.']);
else,
    disp(['Loading file ' datafile ' using NeuroShare library v',...
        num2str(libinfo.LibVersionMaj),'.',num2str(libinfo.LibVersionMin),' ...']);
end;

% open file
[st,fhandle] = ns_OpenFile(datafile);
if st,
    [st,mesg] = ns_GetLastErrorMsg;
    error(mesg);
end;

try,
    % file header
    [st,fheader]  = ns_GetFileInfo(fhandle);
    if st,
        [st,mesg] = ns_GetLastErrorMsg;
        ns_CloseFile(fhandle);
        error(mesg);
    end;

    % Build catalogue of entities
    [st, entityinfo] = ns_GetEntityInfo(fhandle, [1:fheader.EntityCount]);

    % get the channel numbers of analogue and event channels
    % Neuroshare entity types:
    %         Unknown entity                   0
    %         Event entity                     1
    %         Analog entity                    2
    %         Segment entity                   3
    %         Neural event entity              4
    eventlist  = find([entityinfo.EntityType] == 1);
    analoglist = find([entityinfo.EntityType] == 2);

    % How many of a particular entity do we have
    n_event  = length(eventlist);
    n_analog = length(analoglist);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HEADER
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fieldtrip uses:
    %   hdr.Fs           sampling frequency
    %   hdr.nChans       number of channels
    %   hdr.nSamples     number of samples per trial
    %   hdr.nSamplesPre  number of pre-trigger samples in each trial
    %   hdr.nTrials      number of trials
    %   hdr.label        cell-array with labels of each channel
    cnt = 0;
    for channr = analoglist(:)'
        cnt                            = cnt+1;
        [st,einfo] = ns_GetEntityInfo(fhandle,channr);
        [st,ainfo] = ns_GetAnalogInfo(fhandle,channr);
        % Until now, I haven't found a way to check whether the SON files
        % are recorded in continuous or triggered mode by using header
        % information only. Therefore, mode is explicitly set to continuous now.
        % When reading a segment of analog data however, the NeuroShare
        % routines provide information about the continuity of the 
        % read segment. If a discontinuity is detected a warning will
        % be given in the command window.
        if strcmpi(MODE,'continuous')
            out.header(cnt).label        = einfo.EntityLabel;
            out.header(cnt).nsamples     = einfo.ItemCount;
            out.header(cnt).index        = cnt;
            out.header(cnt).sonentityid  = channr; % this is the actual channr in the SON file 
            out.header(cnt).ntrials      = 1;
            out.header(cnt).samplerate   = ainfo.SampleRate;
            out.header(cnt).units        = ainfo.Units;
            out.header(cnt).mode         = 'continuous';
        elseif strcmpi(MODE,'triggered'),
            warning(['Triggered channel mode not implemented yet']);
            out = [];
            return
        else,
            error(['Unknown channel mode for channel ',num2str(channr)]);
        end;
    end;
    out.header = orderfields(out.header);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EVENTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fieldtrip uses:
    %   event.type     = string
    %   event.sample   = expressed in samples, first sample of file is 1
    %   event.value    = number or string
    %   event.offset   = expressed in samples
    %   event.duration = expressed in samples
    if strcmp(pars.readevents,'yes')
        cnt = 0;
        for channr = eventlist(:)'
            [st,einfo]          = ns_GetEntityInfo(fhandle,channr);
            [st,evtime,evdata]  = ns_GetEventData(fhandle,channr,[1:einfo.ItemCount]);

            for trignr = 1:einfo.ItemCount,
                cnt                       = cnt+1;
                out.events(cnt).type      = einfo.EntityLabel;
                out.events(cnt).timestamp = evtime(trignr);
                % Convert timestamp to sample nr. First analog channel's 
                % sample numbers are used as the reference for all
                % timestamps.
                [st,tmp]                  = ns_GetIndexByTime(fhandle,analoglist(1),evtime(trignr),0);
                out.events(cnt).sample    = tmp;
                if iscell(evdata(trignr))
                    out.events(cnt).value = evdata{trignr};  % cell2str
                else
                    out.events(cnt).value = evdata(trignr);
                end;
                out.events(cnt).offset    = 0;
                out.events(cnt).duration  = 0;
            end;
        end;
        out.events = orderfields(out.events);
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if strcmp(pars.readdata,'yes')

          % if no channels specified default to all analog channels
          if isempty(pars.channels), 
              pars.channels = analoglist; 
          else, % renumber requested channels to entityIDs
              if length(pars.channels)>length(analoglist),
                  error(['Requested more analog channels than present in datafile']);
              else,
                  pars.channels = analoglist(pars.channels);
              end;
          end;
                    
          % use first analog channel as reference channel
          targetchan = analoglist(1);           
          [st,einfo] = ns_GetEntityInfo(fhandle,targetchan);
          
          %if no begsample specified default to start of channel
          if pars.begsample<0, begsample=1; else begsample=pars.begsample; end;
          %if no endsample specified default to end of channel
          if pars.endsample<0, endsample=einfo.ItemCount; else endsample=pars.endsample; end;
          
          % calculate number of samples needed (for all requested channels)
          itemcount         = endsample-begsample+1;                  
          targetstep        = out.header([out.header.sonentityid]==targetchan).samplerate;
          [st,targetbegin]  = ns_GetTimeByIndex(fhandle,targetchan,begsample);

          cnt = 0;          
          for channr = pars.channels(:)',
              cnt = cnt+1;
              
              % Cut out the data of channels 1-N in the following way
              %
              %            begtime                        itemcount
              % analog 1 [x..|x..x..x..x..x..x..x..x..x..x..x|..]
              % analog 2 [.x.|.x..x..x..x..x..x..x..x..x..x..x|.]
              % analog N [..x|..x..x..x..x..x..x..x..x..x..x..x|]
              %
              % All channels will have the same number of samples,
              % but timestamps of the channels are offset on subsample
              % level, due to asynchronous sampling of the CED DAQ
              % hardware.

              % return the sample nr occuring after and inclusive of begtime.
              [st,begsample]          = ns_GetIndexByTime(fhandle,channr,targetbegin,1);
              % get the same number of samples as in the target channel
              endsample               = begsample + itemcount-1;
              [st,contcount,chandata] = ns_GetAnalogData(fhandle,channr, begsample,itemcount);
              if contcount~=itemcount, warning(['Discontinuity in data']); end;

              % make a time scale
              [st,begtime]            = ns_GetTimeByIndex(fhandle,channr,begsample);
              [st,endtime]            = ns_GetTimeByIndex(fhandle,channr,endsample);          
              sourcestep              = out.header([out.header.sonentityid]==channr).samplerate;
    
              if sourcestep~=targetstep, warning(['Source and target channels have time steps of different size']); end;
              chantime                = [begtime:1/sourcestep:endtime];
             
              out.data{cnt} = chandata(:)';              
              if strcmp(pars.readtimestamps,'yes') out.time{cnt} = chantime(:)'; end;
          end;
      end;

      % close file
      [st]  = ns_CloseFile(fhandle);
    if st,
        [st,mesg] = ns_GetLastErrorMsg;
        disp(mesg);
    end;
   
% use catch to close any opened files before terminating
catch,
    %   fclose(fid);
    [st]  = ns_CloseFile(fhandle);
    if st,
        [st,mesg] = ns_GetLastErrorMsg;
        disp(mesg);
    end;
    error(lasterr);
end;

