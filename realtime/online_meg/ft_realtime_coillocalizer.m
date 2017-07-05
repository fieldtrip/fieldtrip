function ft_realtime_coillocalizer(cfg)

% FT_REALTIME_COILLOCALIZER is a realtime application for online tracking
% of MEG localizer coils.
%
% Use as
%   ft_realtime_coillocalizer(cfg)
% with the following configuration options
%   cfg.blocksize  = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.channel    = cell-array, see FT_CHANNELSELECTION (default = {'MEG', 'MEGREF'})
%   cfg.bufferdata = whether to process the 'first or 'last' data that is available (default = 'last')
%   cfg.jumptoeof  = whether to skip to the end of the stream/file at startup (default = 'yes')
%
% The settings for extracting the spatial topgraphy of each coil are configured as
%   cfg.coilfreq      = single number in Hz or list of numbers
%   cfg.refchan       = single string or cell-array with strings
%
% The source of the data is configured as
%   cfg.dataset       = string
% or alternatively to obtain more low-level control as
%   cfg.datafile      = string
%   cfg.headerfile    = string
%   cfg.eventfile     = string
%   cfg.dataformat    = string, default is determined automatic
%   cfg.headerformat  = string, default is determined automatic
%   cfg.eventformat   = string, default is determined automatic
%
% Some notes about skipping data and catching up with the data stream:
%
% cfg.jumptoeof='yes' causes the realtime function to jump to the end
% when the function _starts_. It causes all data acquired prior to
% starting the RT function to be skipped.
%
% cfg.bufferdata=last causes the realtime function to jump to the last
% available data while _running_. If the realtime loop is not fast enough,
% it causes some data to be dropped.
%
% If you want to skip all data that was acquired before you start the
% realtime function, but don't want to miss any data that was acquired while
% the realtime function is started, then you should use jumptoeof=yes and
% bufferdata=first. If you want to analyse data from a file, then you
% should use jumptoeof=no and bufferdata=first.
%
% To stop this realtime function, you have to press Ctrl-C

% Copyright (C) 2011-2012, Robert Oostenveld
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

% set the default configuration options
if ~isfield(cfg, 'dataformat'),     cfg.dataformat = [];      end % default is detected automatically
if ~isfield(cfg, 'headerformat'),   cfg.headerformat = [];    end % default is detected automatically
if ~isfield(cfg, 'eventformat'),    cfg.eventformat = [];     end % default is detected automatically
if ~isfield(cfg, 'blocksize'),      cfg.blocksize = 1;        end % in seconds
if ~isfield(cfg, 'overlap'),        cfg.overlap = 0;          end % in seconds
if ~isfield(cfg, 'channel'),        cfg.channel = {'MEG', 'MEGREF'}; end
if ~isfield(cfg, 'refchan'),        cfg.refchan = {};         end
if ~isfield(cfg, 'coilfreq'),       cfg.coilfreq = [];        end
if ~isfield(cfg, 'bufferdata'),     cfg.bufferdata = 'last';  end % first or last
if ~isfield(cfg, 'jumptoeof'),      cfg.jumptoeof = 'no';     end % jump to end of file at initialization

if ~iscell(cfg.refchan)
    % convert from string to cell-array
    cfg.refchan = {cfg.refchan};
end

if ~isfield(cfg, 'dataset') && ~isfield(cfg, 'header') && ~isfield(cfg, 'datafile')
    cfg.dataset = 'buffer://localhost:1972';
end

% translate dataset into datafile+headerfile
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');
cfg = ft_checkconfig(cfg, 'required', {'datafile' 'headerfile'});

% ensure that the persistent variables related to caching are cleared
clear ft_read_header

% start by reading the header from the realtime buffer
hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true, 'retry', true);

isneuromag = ft_senstype(hdr.grad, 'neuromag306');
isctf      = ft_senstype(hdr.grad, 'ctf275');

if isneuromag
    ctr = 0;
    for j = 1:length(hdr.orig.raw.info.dig)
        if hdr.orig.raw.info.dig(1,j).kind == 2 % 1: RPA/LPA/NAS, 2: HPI coils, 4: Head shape
            ctr = ctr+1;
            hc(ctr,:) = double((hdr.orig.raw.info.dig(1,j).r).*100); % convert from m to cm
        end
    end
    fprintf('%d HPI coils found in fif header\n', ctr);
end


% define a subset of channels for reading
cfg.channel = ft_channelselection(cfg.channel, hdr.label);
cfg.refchan = ft_channelselection(cfg.refchan, hdr.label);
megindx     = match_str(hdr.label, cfg.channel);
refindx     = match_str(hdr.label, cfg.refchan);
chanindx    = sort([megindx(:)' refindx(:)']);

% sofar the megindx and refindx were indices into all data channels
% change them into indexing vectors into the selected channels
[dum, megindx] = intersect(chanindx, megindx);
[dum, refindx] = intersect(chanindx, refindx);

nchan = length(chanindx);
if nchan==0
    ft_error('no channels were selected');
end

% determine the size of blocks to process
blocksize = round(cfg.blocksize * hdr.Fs);
overlap   = round(cfg.overlap*hdr.Fs);

% construct the reference signal for each of the coils
% FIXME the blocksize should match an integer number of cycles -> perhaps use nan and nansum?
ncoil = length(cfg.coilfreq);
if ncoil==0
    ft_error('no coil frequencies were specified');
else
    time = (1:blocksize)./hdr.Fs;
    coil = zeros(ncoil, blocksize);
    for i=1:ncoil
        coil(i,:) = exp(time*cfg.coilfreq(i)*1i*2*pi);
        coil(i,:) = coil(i,:) / norm(coil(i,:));
    end
end

% prepare the forward model and the sensor array for subsequent fitting
% note that the forward model is a magnetic dipole in an infinite vacuum
[vol, sens] = ft_prepare_vol_sens([], hdr.grad, 'channel', cfg.channel);

% open a figure
figure

% set an initial guess for each of the dipole/coil positions
if isneuromag
    for i=1:ncoil
        dip(i).pos = hc(i,:);
        dip(i).mom = [0 0 0]';
    end
else
    for i=1:ncoil
        dip(i).pos = mean(sens.coilpos,1); % somewhere in the middle of the helmet
        dip(i).mom = [0 0 0]';
    end
end

if strcmp(cfg.jumptoeof, 'yes')
    prevSample = hdr.nSamples * hdr.nTrials;
elseif isfield(cfg, 'offset')
    prevSample = (cfg.offset*hdr.Fs);
else
    prevSample  = 0;
end
count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while true
    
    % determine number of samples available in buffer
    hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat, 'cache', true);
    
    % see whether new samples are available
    newsamples = (hdr.nSamples*hdr.nTrials-prevSample);
    
    if newsamples>=blocksize
        
        % determine the samples to process
        if strcmp(cfg.bufferdata, 'last')
            begsample  = hdr.nSamples*hdr.nTrials - blocksize + 1;
            endsample  = hdr.nSamples*hdr.nTrials;
        elseif strcmp(cfg.bufferdata, 'first')
            begsample  = prevSample+1;
            endsample  = prevSample+blocksize ;
        else
            ft_error('unsupported value for cfg.bufferdata');
        end
        
        % this allows overlapping data segments
        if overlap && (begsample>overlap)
            begsample = begsample - overlap;
            endsample = endsample - overlap;
        end
        
        % remember up to where the data was read
        prevSample  = endsample;
        count       = count + 1;
        fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);
        
        % read data segment from buffer
        dat = ft_read_data(cfg.datafile, 'header', hdr, 'dataformat', cfg.dataformat, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % from here onward it is specific to the processing of the data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % split the data into the MEG channels and the (optional) reference channels
        meg = dat(megindx,:);
        ref = dat(refindx,:);
        
        if ~isempty(ref)
            % only consider the signal that is phase-locked to the reference channels
            ampl  = sum(ref.*coil);   % estimate the complex amplitude
            phase = ampl./abs(ampl);  % estimate the phase
            % FIXME the code should be checked for phase-flips
            for i=1:ncoil
                topo(:,i) = topo(:,i) * phase(i);
            end
        else
            % estimate the complex-valued MEG topography for each coil
            % this implements a discrete Fourier transform (DFT)
            topo = ft_preproc_detrend(meg) * ctranspose(coil);
        end
        
        % ignore the out-of-phase spectral component in the topography
        topo = real(topo); % THIS SEEMS TO BE CRUCIAL
        
        if false
            close all
            for i=1:ncoil
                figure
                ft_plot_sens(sens);
                ft_plot_dipole(dip(i).pos, dip(i).mom);
                ft_plot_topo3d(sens.chanpos, real(topo(:,i)));
                m = max(abs(caxis));
                caxis([-m m]);
                alpha 0.5
            end
        end
        
        % fit a magnetic dipole to each of the topographies
        if isneuromag
            constr.sequential = true;
            constr.rigidbody = true;
            % fit the coils together
            for i=1:ncoil
                pos(i,:) = dip(i).pos;
            end
            dipall = [];
            dipall.pos = pos;
            dipall = dipole_fit(dipall, sens, vol, topo, 'constr', constr);
            for i=1:ncoil
                sel = (1:3) + 3*(i-1);
                dip(i).pos = dipall.pos(i,:);
                dip(i).mom = real(dipall.mom(sel,i)); % ignore the complex phase information
            end
        else
            % fit the coils sequentially
            for i=1:ncoil
                dip(i) = dipole_fit(dip(i), sens, vol, topo(:,i));
            end
        end
        
        cla
        % plot the gradiometer array for reference
        ft_plot_sens(sens);
        % plot each of the fitted dipoles
        for i=1:ncoil
            ft_plot_dipole(dip(i).pos, dip(i).mom);
        end
        
        % show current timesample
        str = sprintf('samples: %d - %d\n', begsample, endsample);
        title(str);
        
        % force Matlab to update the figure
        drawnow
        
    end % if enough new samples
end % while true
