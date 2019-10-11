% eeglab2fieldtrip() - This function converts EEGLAB datasets to Fieldtrip
%                      for source localization (DIPFIT).
%
% Usage:
%  >> data = eeglab2fieldtrip( EEG, fieldbox, transform );
%
% Inputs:
%   EEG       - [struct] EEGLAB structure
%   fieldbox  - ['raw'|'timelock'|'comp'] for reasons, you can also
%               specify 'preprocessing' (raw), 'timelockanalysis' (timelock)
%               and 'componentanalysis' (comp).
%   transform - ['none'|'dipfit'] transform channel locations for DIPFIT
%               using the transformation matrix in the field
%               'coord_transform' of the dipfit substructure of the EEG
%               structure.
% Outputs:
%   data      - FIELDTRIP data structure (see data type functions below)
%
% Note:
%   For conversion of low level raw data structures, save the EEGLAB
%   dataset and use:
%      hdr    = ft_read_header( EEGLABFILE ); 
%      data   = ft_read_data(   EEGLABFILE, 'header', hdr ); 
%      events = ft_read_event(  EEGLABFILE, 'header', hdr );
%
% Author: Robert Oostenveld, F.C. Donders Centre, May, 2004.
%         Arnaud Delorme, SCCN, INC, UCSD
%
% See also: FT_DATATYPE_RAW, FT_DATATYPE_TIMELOCK, FT_DATATYPE_COMP

% Copyright (C) 2004 Robert Oostenveld, F.C. Donders Centre, roberto@smi.auc.dk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function data = eeglab2fieldtrip(EEG, fieldbox, transform)

if nargin < 2
    help eeglab2fieldtrip
    return;
end

% start with an empty data object
data = [];

% add the objects that are common to all fieldboxes
tmpchanlocs  = EEG.chanlocs;

% get the electrode positions from the EEG structure: in principle, the number of
% channels can be more or less than the number of channel locations, i.e. not
% every channel has a position, or the potential was not measured on every
% position. This is not supported by EEGLAB, but it is supported by FIELDTRIP.

if strcmpi(fieldbox, 'chanloc_withfid')
    % insert "no data channels" in channel structure
    % ----------------------------------------------
    if isfield(EEG.chaninfo, 'nodatchans') && ~isempty( EEG.chaninfo.nodatchans )
        chanlen = length(EEG.chanlocs);
        fields = fieldnames( EEG.chaninfo.nodatchans );
        for index = 1:length(EEG.chaninfo.nodatchans)
            ind = chanlen+index;
            for f = 1:length( fields )
                EEG.chanlocs = setfield(EEG.chanlocs, { ind }, fields{f}, ...
                    getfield( EEG.chaninfo.nodatchans, { index },  fields{f}));
            end
        end
    end
end

if ~isempty(EEG.chanlocs) && ~isempty(EEG.chanlocs(1).X)
    data.elec.elecpos = zeros(length( EEG.chanlocs ), 3);
    for ind = 1:length( EEG.chanlocs )
        data.elec.label{ind} = EEG.chanlocs(ind).labels;
        if ~isempty(EEG.chanlocs(ind).X)
            data.elec.elecpos(ind,1) = EEG.chanlocs(ind).X;
            data.elec.elecpos(ind,2) = EEG.chanlocs(ind).Y;
            data.elec.elecpos(ind,3) = EEG.chanlocs(ind).Z;
        else
            data.elec.elecpos(ind,:) = [0 0 0];
        end
    end
    data.elec.pnt = data.elec.elecpos;
end

if nargin > 2
    if strcmpi(transform, 'dipfit')
        if ~isempty(EEG.dipfit.coord_transform)
            disp('Transforming electrode coordinates to match head model');
            transfmat = traditionaldipfit(EEG.dipfit.coord_transform);
            data.elec.elecpos = transfmat * [ data.elec.elecpos ones(size(data.elec.elecpos,1),1) ]';
            data.elec.elecpos = data.elec.elecpos(1:3,:)';
            data.elec.pnt = data.elec.elecpos;
        else
            disp('Warning: no transformation of electrode coordinates to match head model');
        end
    end
end

switch fieldbox
    case { 'preprocessing' 'raw' }
        data.fsample = EEG.srate;
        for index = 1:EEG.trials
            data.trial{index}  = EEG.data(:,:,index);
            data.time{index}   = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
        end
        data.label = getchanlabels(tmpchanlocs, 1:EEG.nbchan);
        if EEG.trials > 1
            res = std_maketrialinfo([], EEG);
            data.trialinfo = struct2table(res.datasetinfo.trialinfo);
        else
            res = struct('subject', EEG.subject, 'condition', EEG.condition, 'group', EEG.group, 'session', EEG.session);
            data.trialinfo = struct2table(res);
            if isempty(data.trialinfo)
                data = rmfield(data, 'trialinfo');
            end
        end
        
    case { 'timelock' 'timelockanalysis' }
        data.avg  = mean(EEG.data, 3);
        data.var  = std(EEG.data, [], 3).^2;
        data.time = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
        data.label = getchanlabels(tmpchanlocs, 1:EEG.nbchan);
        
    case { 'comp' 'componentanalysis' }
        if isempty(EEG.icaact)
            icaacttmp = eeg_getica(EEG);
        end
        for index = 1:EEG.trials
            % the trials correspond to the raw data trials, except that they
            % contain the component activations
            try
                if isempty(EEG.icaact)
                    data.trial{index} = icaacttmp(:,:,index); % Using icaacttmp to not change EEG structure
                else
                    data.trial{index}  = EEG.icaact(:,:,index);
                end
            catch
                
            end
            data.time{index}   = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
        end
        data.label = [];
        for comp = 1:size(EEG.icawinv,2)
            % the labels correspond to the component activations that are stored in data.trial
            data.label{comp} = sprintf('ica_%03d', comp);
        end
        % get the spatial distribution and electrode positions
        tmpchanlocs    = EEG.chanlocs;
        data.topolabel = { tmpchanlocs(EEG.icachansind).labels };
        data.topo      = EEG.icawinv;
		% Copy weights & sphere, too
    	data.unmixing = EEG.icaweights*EEG.icasphere;
        
    case { 'chanloc' 'chanloc_withfid' }
        
    otherwise
        error('unsupported fieldbox')
end

try
    % get the full name of the function
    data.cfg.version.name = mfilename('fullpath');
catch
    % required for compatibility with Matlab versions prior to release 13 (6.5)
    [st, i] = dbstack;
    data.cfg.version.name = st(i);
end

% convert EEGLAB channel labels to Fieldtrip channel labels
function label = getchanlabels(tmpchanlocs, indices)
if ~isempty(tmpchanlocs)
    label   = { tmpchanlocs.labels };
else
    for iChan = indices
        label{iChan} = num2str(iChan);
    end
end


