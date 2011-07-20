function ft_neighbourplot(cfg, data)

% FT_NEIGHBOURPLOT visualizes neighbouring channels in a particular channel
% configuration. The positions of the channel are specified in a
% gradiometer or electrode configuration or from a layout. Neighbouring
% channels are obtained by ft_neighbourselection.
%
% Use as
%  ft_neighbourplot(cfg)
% or as
%  ft_neighbourplot(cfg, data)
%
% where
%
%   cfg.neighbours    = neighbour structure from ft_neighbourselection
%   (optional)
%   cfg.elec          = structure with EEG electrode positions
%   cfg.grad          = structure with MEG gradiometer positions
%   cfg.elecfile      = filename containing EEG electrode positions
%   cfg.gradfile      = filename containing MEG gradiometer positions
%   cfg.layout        = filename of the layout, see FT_PREPARE_LAYOUT
%   cfg.verbose       = 'yes' or 'no', if 'yes' then the plot callback will
%                       include text output
%
% The following data fields may also be used by FT_NEIGHBOURSELECTION:
%   data.elec     = structure with EEG electrode positions
%   data.grad     = structure with MEG gradiometer positions
%
%
% If cfg.neighbours is no defined or empty, the function calls 
% ft_neighbourselection to compute channel neighbours. 
%
% See also FT_NEIGHBOURSELECTION


% Copyright (C) 2011, Jörn M. Horschig, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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


if nargin < 3
    if isfield(cfg, 'neighbours')
        neighbours = cfg.neighbours;
    elseif nargin < 2
        neighbours = ft_neighbourselection(cfg);
    else
        neighbours = ft_neighbourselection(cfg, data);
    end
end


if iscell(cfg.neighbours)
    warning('Neighbourstructure is in old format - converting to structure array');
    cfg.neighbours = fixneighbours(cfg.neighbours);
end

if ~isfield(cfg, 'verbose')
    cfg.verbose = 'no';
elseif strcmp(cfg.verbose, 'yes')    
    cfg.verbose = true;
end

% get the the grad or elec if not present in the data
if isfield(data, 'grad')
    fprintf('Using the gradiometer configuration from the dataset.\n');
    sens = data.grad;
    % extract true channelposition
    [sens.pnt, sens.ori, sens.label] = channelposition(sens);
elseif isfield(data, 'elec')
    fprintf('Using the electrode configuration from the dataset.\n');
    sens = data.elec;
elseif isfield(cfg, 'grad')
    fprintf('Obtaining the gradiometer configuration from the configuration.\n');
    sens = cfg.grad;
    % extract true channelposition
    [sens.pnt, sens.ori, sens.label] = channelposition(sens);
elseif isfield(cfg, 'elec')
    fprintf('Obtaining the electrode configuration from the configuration.\n');
    sens = cfg.elec;
elseif isfield(cfg, 'gradfile')
    fprintf('Obtaining the gradiometer configuration from a file.\n');
    sens = ft_read_sens(cfg.gradfile);
    % extract true channelposition
    [sens.pnt, sens.ori, sens.label] = channelposition(sens);
elseif isfield(cfg, 'elecfile')
    fprintf('Obtaining the electrode configuration from a file.\n');
    sens = ft_read_sens(cfg.elecfile);
elseif isfield(cfg, 'layout')
    fprintf('Using the 2-D layout to determine the neighbours\n');
    lay = ft_prepare_layout(cfg);
    sens = [];
    sens.label = lay.label;
    sens.pnt = lay.pos;
    sens.pnt(:,3) = 0;
else
    error('Did not find gradiometer or electrode information.');
end;

% give some graphical feedback
if all(sens.pnt(:,3)==0)
    % the sensor positions are already projected on a 2D plane
    proj = sens.pnt(:,1:2);
else
    % use 3-dimensional data for plotting
    proj = sens.pnt;
end
figure
axis equal
axis off
hold on;
for i=1:length(neighbours)
    this = neighbours(i);    
    
    sel1 = match_str(sens.label, this.label);
    sel2 = match_str(sens.label, this.neighblabel);
    
    for j=1:length(this.neighblabel)
        x1 = proj(sel1,1);
        y1 = proj(sel1,2);
        x2 = proj(sel2(j),1);
        y2 = proj(sel2(j),2);
        X = [x1 x2];
        Y = [y1 y2];
        if size(proj, 2) == 2
            line(X, Y, 'color', 'r');
        elseif size(proj, 2) == 3
            z1 = proj(sel1,3);
            z2 = proj(sel2(j),3);
            Z = [z1 z2];
            line(X, Y, Z, 'color', 'r');
        end
    end
end

% this is for putting the channels on top of the connections
set(gcf, 'UserData', []);
for i=1:length(neighbours)
    if size(proj, 2) == 2
        hs(i) = line(proj(i, 1), proj(i, 2),                                            ...     
                        'MarkerEdgeColor',  'k',                                        ...
                        'MarkerFaceColor',  'k',                                        ...
                        'Marker',           'o',                                        ...
                        'MarkerSize',       .125*(2+numel(neighbours(i).neighblabel)^2), ...
                        'UserData',         i,                                          ...
                        'ButtonDownFcn',    @showLabelInTitle);
                    
    elseif size(proj, 2) == 3
        hs(i) = line(proj(i, 1), proj(i, 2), proj(i, 3),                                ...
                        'MarkerEdgeColor',  'k',                                        ...
                        'MarkerFaceColor',  'k',                                        ...
                        'Marker',           'o',                                        ...
                        'MarkerSize',       .125*(2+numel(neighbours(i).neighblabel)^2), ...
                        'UserData',         i,                        ...
                        'ButtonDownFcn',    @showLabelInTitle);
    else
        error('Channel coordinates are too high dimensional');
    end
end
hold off;
title('[Click on a sensor to see its label]');

    function showLabelInTitle(gcbo, EventData, handles)
        lastSensId  = get(gcf,  'UserData');
        curSensId   = get(gcbo, 'UserData');
        
        title(['Selected channel: ' neighbours{curSensId}.label]);
        if cfg.verbose 
            fprintf('Selected channel: %s\n', neighbours{curSensId}.label); 
        end
        
        set(hs(curSensId), 'MarkerFaceColor', 'g');
        set(hs(lastSensId), 'MarkerFaceColor', 'k');
        
        set(gcf, 'UserData', curSensId');
    end
end