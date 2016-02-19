function ft_realtime_fileproxy(cfg)

% FT_REALTIME_FILEPROXY reads continuous data from an EEG/MEG file and writes it to a
% FieldTrip buffer. This works for any file format that is supported by FieldTrip.
%
% The FieldTrip buffer is a network transparent server that allows the acquisition
% client to stream data to it. An analysis client can connect to read the data upon
% request. Multiple clients can connect simultaneously, each analyzing a specific
% aspect of the data concurrently.
%
% Use as
%   ft_realtime_fileproxy(cfg)
% with the following configuration options
%   cfg.minblocksize         = number, in seconds (default = 0)
%   cfg.maxblocksize         = number, in seconds (default = 1)
%   cfg.channel              = cell-array, see FT_CHANNELSELECTION (default = 'all')
%   cfg.jumptoeof            = jump to end of file at initialization (default = 'no')
%   cfg.readevent            = whether or not to copy events (default = 'no'; event type can also be specified; e.g., 'UPPT002')
%   cfg.speed                = relative speed at which data is written (default = inf)
%
% The source of the data is configured as
%   cfg.source.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.source.datafile      = string
%   cfg.source.headerfile    = string
%   cfg.source.eventfile     = string
%   cfg.source.dataformat    = string, default is determined automatic
%   cfg.source.headerformat  = string, default is determined automatic
%   cfg.source.eventformat   = string, default is determined automatic
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C
%
% See also FT_REALTIME_SIGNALPROXY, FT_REALTIME_SIGNALVIEWER

% Copyright (C) 2008, Robert Oostenveld
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

% set the defaults
if ~isfield(cfg, 'source'),               cfg.source = [];                                  end
if ~isfield(cfg, 'target'),               cfg.target = [];                                  end
if ~isfield(cfg.source, 'headerformat'),  cfg.source.headerformat = [];                     end % default is detected automatically
if ~isfield(cfg.source, 'dataformat'),    cfg.source.dataformat = [];                       end % default is detected automatically
if ~isfield(cfg.target, 'headerformat'),  cfg.target.headerformat = [];                     end % default is detected automatically
if ~isfield(cfg.target, 'dataformat'),    cfg.target.dataformat = [];                       end % default is detected automatically
if ~isfield(cfg.target, 'datafile'),      cfg.target.datafile = 'buffer://localhost:1972';  end
if ~isfield(cfg, 'minblocksize'),         cfg.minblocksize = 0;                             end % in seconds
if ~isfield(cfg, 'maxblocksize'),         cfg.maxblocksize = 1;                             end % in seconds
if ~isfield(cfg, 'channel'),              cfg.channel = 'all';                              end
if ~isfield(cfg, 'jumptoeof'),            cfg.jumptoeof = 'no';                             end % jump to end of file at initialization
if ~isfield(cfg, 'readevent'),            cfg.readevent = 'no';                             end % capture events?
if ~isfield(cfg, 'speed'),                cfg.speed = inf ;                                 end % inf -> run as fast as possible

% translate dataset into datafile+headerfile
cfg.source = ft_checkconfig(cfg.source, 'dataset2files', 'yes');
cfg.target = ft_checkconfig(cfg.target, 'dataset2files', 'yes');
ft_checkconfig(cfg.source, 'required', {'datafile' 'headerfile'});
ft_checkconfig(cfg.target, 'required', {'datafile' 'headerfile'});

if ~isfield(cfg.source,'eventfile') || isempty(cfg.source.eventfile)
    cfg.source.eventfile = cfg.source.datafile;
end

if ~isfield(cfg.target,'eventfile') || isempty(cfg.target.eventfile)
    cfg.target.eventfile = cfg.target.datafile;
end

% ensure that the persistent variables related to caching are cleared
clear ft_read_header
% read the header for the first time
hdr = ft_read_header(cfg.source.headerfile);
fprintf('updating the header information, %d samples available\n', hdr.nSamples*hdr.nTrials);

% define a subset of channels for reading
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
chanindx    = match_str(hdr.label, cfg.channel);
nchan       = length(chanindx);
if nchan==0
    error('no channels were selected');
end

minblocksmp = round(cfg.minblocksize*hdr.Fs);
minblocksmp = max(minblocksmp, 1);
maxblocksmp = round(cfg.maxblocksize*hdr.Fs);
count       = 0;

if strcmp(cfg.jumptoeof, 'yes')
    prevSample = hdr.nSamples * hdr.nTrials;
else
    prevSample  = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evt = [];
while true

    % determine number of samples available in buffer
    hdr = ft_read_header(cfg.source.headerfile, 'cache', true);

    % see whether new samples are available
    newsamples = (hdr.nSamples*hdr.nTrials-prevSample);

    if newsamples>=minblocksmp

        begsample  = prevSample+1;
        endsample  = prevSample+min(newsamples,maxblocksmp);

        % remember up to where the data was read
        count       = count + 1;
        fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);

        % read data segment
        dat = ft_read_data(cfg.source.datafile,'header', hdr, 'dataformat', cfg.source.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);

        % it only makes sense to read those events associated with the currently processed data
        if ~strcmp(cfg.readevent,'no')
          evt = ft_read_event(cfg.source.eventfile, 'header', hdr, 'minsample', begsample, 'maxsample', endsample);

          if ~strcmp(cfg.readevent,'yes')
            evt = ft_filter_event(evt, 'type', cfg.readevent);
          end
          
        end

        prevSample  = endsample;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % from here onward it is specific to writing the data to another stream
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if count==1
            % the input file may have a different offset than the output file
            offset = begsample - 1;
            % flush the file, write the header and subsequently write the data segment
            ft_write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', false);
            if ~strcmp(cfg.readevent,'no')
                for i=1:numel(evt)
                    evt(i).sample = evt(i).sample - offset;
                end
                ft_write_event(cfg.target.eventfile,evt,'append',false);
            end
        else
            % write the data segment
            ft_write_data(cfg.target.datafile, dat, 'header', hdr, 'dataformat', cfg.target.dataformat, 'chanindx', chanindx, 'append', true);
            if ~strcmp(cfg.readevent,'no')
                for i=1:numel(evt)
                    evt(i).sample = evt(i).sample - offset;
                end
                ft_write_event(cfg.target.eventfile,evt,'append',true);
            end
        end % if count==1

        % wait for a realistic amount of time
        pause(((endsample-begsample+1)/hdr.Fs)/cfg.speed);

    end % if enough new samples
end % while true
