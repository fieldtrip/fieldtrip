% eeglab2fieldtrip() - This function converts EEGLAB datasets to Fieldtrip
%                      for source localization (DIPFIT). IN PRACTICE TO 
%                      CONVERT DATASETS, IT IS RECOMMENDED TO USE FILE-IO
%                      see https://sccn.ucsd.edu/wiki/EEGLAB_and_Fieldtrip
%
% Usage:
%  >> [header,dat,evt] = eeglab2fieldtrip( EEG, fieldbox, transform );
%
% Inputs:
%   EEG       - [struct] EEGLAB structure
%   fieldbox  - ['preprocessing'|'freqanalysis'|'timelockanalysis'|'componentanalysis']
%   transform - ['none'|'dipfit'] transform channel locations for DIPFIT
%               using the transformation matrix in the field
%               'coord_transform' of the dipfit substructure of the EEG
%               structure.
% Outputs:
%   header    - FIELDTRIP header structure
%   dat       - raw data (compatible with Fieldtrip)
%   evt       - FIELDTRIP event structure
%
% Author: Robert Oostenveld, F.C. Donders Centre, May, 2004.
%         Arnaud Delorme, SCCN, INC, UCSD
%
% See also:

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

function [header,data,evt] = eeglab2fieldtrip(EEG, fieldbox, transform)

if nargin < 2
    help eeglab2fieldtrip
    return;
end;

% start with an empty header object
header = [];

% add the objects that are common to all fieldboxes
tmpchanlocs  = EEG.chanlocs;
header.label   = { tmpchanlocs(EEG.icachansind).labels };
header.fsample = EEG.srate;

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
            end;
        end;
    end;
end;

header.elec.pnt   = zeros(length( EEG.chanlocs ), 3);
for ind = 1:length( EEG.chanlocs )
    header.elec.label{ind} = EEG.chanlocs(ind).labels;
    if ~isempty(EEG.chanlocs(ind).X)
        header.elec.pnt(ind,1) = EEG.chanlocs(ind).X;
        header.elec.pnt(ind,2) = EEG.chanlocs(ind).Y;
        header.elec.pnt(ind,3) = EEG.chanlocs(ind).Z;
    else
        header.elec.pnt(ind,:) = [0 0 0];
    end;
end;

if nargin > 2
    if strcmpi(transform, 'dipfit')
        if ~isempty(EEG.dipfit.coord_transform)
            disp('Transforming electrode coordinates to match head model');
            transfmat = traditionaldipfit(EEG.dipfit.coord_transform);
            header.elec.pnt = transfmat * [ header.elec.pnt ones(size(header.elec.pnt,1),1) ]';
            header.elec.pnt = header.elec.pnt(1:3,:)';
        else
            disp('Warning: no transformation of electrode coordinates to match head model');
        end;
    end;
end;

switch fieldbox
    case 'preprocessing'
        for index = 1:EEG.trials
            header.trial{index}  = EEG.data(:,:,index);
            header.time{index}   = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
        end;
        header.label   = { tmpchanlocs(1:EEG.nbchan).labels };
        
        
    case 'timelockanalysis'
        header.avg  = mean(EEG.data, 3);
        header.var  = std(EEG.data, [], 3).^2;
        header.time = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
        header.label   = { tmpchanlocs(1:EEG.nbchan).labels };
        
    case 'componentanalysis'
        if isempty(EEG.icaact)
            icaacttmp = eeg_getica(EEG);
        end
        for index = 1:EEG.trials
            % the trials correspond to the raw header trials, except that they
            % contain the component activations
            try
                if isempty(EEG.icaact)
                    header.trial{index} = icaacttmp(:,:,index); % Using icaacttmp to not change EEG structure
                else
                    header.trial{index}  = EEG.icaact(:,:,index);
                end
            catch
                
            end;
            header.time{index}   = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
        end;
        header.label = [];
        for comp = 1:size(EEG.icawinv,2)
            % the labels correspond to the component activations that are stored in header.trial
            header.label{comp} = sprintf('ica_%03d', comp);
        end
        % get the spatial distribution and electrode positions
        tmpchanlocs    = EEG.chanlocs;
        header.topolabel = { tmpchanlocs(EEG.icachansind).labels };
        header.topo      = EEG.icawinv;
		% Copy weights & sphere, too
    	header.weights = EEG.icaweights;
    	header.sphere  = EEG.icasphere;
        
    case { 'chanloc' 'chanloc_withfid' }
        
    case 'freqanalysis'
        error('freqanalysis fieldbox not implemented yet')
        
    otherwise
        error('unsupported fieldbox')
end

try
    % get the full name of the function
    header.cfg.version.name = mfilename('fullpath');
catch
    % required for compatibility with Matlab versions prior to release 13 (6.5)
    [st, i] = dbstack;
    header.cfg.version.name = st(i);
end

data = EEG.data;

% convert event structure

