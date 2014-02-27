function [cfg] = ft_multiplotCC(cfg, data)

% FT_MULTIPLOTCC visualises the connectiivty between channels by using
% multiple topoplots. The topoplot at a given channel location shows the
% coherence of that channel with all other channels.
%
% Use as
%   ft_multiplotCC(cfg, data)
%
% See also FT_PREPARE_LAYOUT, FT_TOPOPLOTCC, FT_CONNECTIVITYPLOT

% Configuration options:
% cfg.layout    = layout filename or a structure produced by prepare_layout
% (requred)
% cfg.foi       = lower and uper lower limits of freqency of interest (if freq is
%                   present) (default = 'all');
% cfg.toi       = lower and uper lower limits of time of interest (if time is
%                   present)   (default = 'all');
% cfg.parameter = the paramater to be plotted. Must be bivariate data.
%                   (dafault = 'cohspctrm')
% cfg.zlim      = plotting limits for color dimension, 'maxmin', 
%                   'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
% cfg.colorbar  = 'off' or 'on', displays colorbar (default = 'off').
% cfg.marker    = 'labels', 'numbers' or 'off' displays channel names next
%                   to each plot (default = 'markers')
% cfg.directionality     = '', 'inflow' or 'outflow' specifies for
%                            connectivity measures whether the inflow into a
%                            node, or the outflow from a node is plotted.
%                            Leave empty for symmetrical (uundirected)
%                            connectivity (default = '');
% cfg.outline   = 'off' or 'on' plots the head outline for each plot
%                   (default = 'off');
%
% Additonal low-level plotting options (see FT_PLOT_TOPO):
% cfg.interpmethod
% cfg.shading
% cfg.gridscale
% cfg.style
% cfg.interplim

% This function requires input from FT_FREQSTATISTICS_SHIFTPREDICT
% This function should be rewritten, using the clean topoplot implementation

% Copyright (C) 2005-2006, Jan-Mathijs Schoffelen, Robert Oostenveld
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
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug

% check if the input data is valid for this function
data = ft_checkdata(data);

% check if the input cfg is valid for this function

if ~isfield(cfg, 'layout'),
error('cfg.layout is required.');
end;
if ~isfield(cfg, 'foilim'),         cfg.foilim   = 'all';               end;
if ~isfield(cfg, 'toilim'),         cfg.toilim   = 'all';               end;
if ~isfield(cfg, 'parameter'),      cfg.parameter = 'cohspctrm';        end;
if ~isfield(cfg, 'zlim'),           cfg.zlim = 'maxmin';                end;
if ~isfield(cfg, 'directionality'), cfg.directionality = [];            end;
if ~isfield(cfg, 'outline'),        cfg.outline = 'off';                 end;

% options for low-level ft_plot_topo function
if ~isfield(cfg, 'interpmethod'),   cfg.interpmethod = 'v4';            end;
if ~isfield(cfg, 'shading'),        cfg.shading = 'flat';               end;
if ~isfield(cfg, 'gridscale'),      cfg.gridscale = 50;                 end;
if ~isfield(cfg, 'style'),          cfg.style = 'surf';                 end;
if ~isfield(cfg, 'interplim'),      cfg.interplim = 'mask';             end;

% options for marking channels
if ~isfield(cfg, 'marker'),         cfg.marker = 'labels';                 end;
if ~isfield(cfg, 'markersymbol'),   cfg.markersymbol = 'o';             end;
if ~isfield(cfg, 'markercolor'),    cfg.markercolor = [0 0 0];          end;
if ~isfield(cfg, 'markersize'),     cfg.markersize = 2;                 end;
if ~isfield(cfg, 'markerfontsize'), cfg.markerfontsize = 8;             end;
if ~isfield(cfg, 'labeloffset'),    cfg.labeloffset = 0.005;            end;



%these arent implimented yet
if ~isfield(cfg, 'channellocal'),   cfg.channellocal = 'all';  end;
if ~isfield(cfg, 'channelglobal'),   cfg.channelglobal = 'all';  end;
if ~isfield(cfg, 'chansubsample'),   cfg.chansubsample = 'off';  end;


% identify dimensions and indices ot average data over
dimtok = tokenize(data.dimord, '_');

idx_string=[];

for i=1:length(dimtok);
    if strcmp(dimtok{i},'chan');
        chan_dim=i;
        idx_string=[idx_string ':,'];
        
    elseif strcmp(dimtok{i},'freq');
        fdim=i;
        if strcmp(cfg.foilim, 'all');
            foilim=[1 length(data.freq)];
        else
            foilim=[nearest(data.freq, cfg.foilim(1)) nearest(data.freq, cfg.foilim(end))];
        end
        
        idx_string=[idx_string num2str(foilim(1)) ':' num2str(foilim(end)) ','];
        
        
        
        
    elseif strcmp(dimtok{i},'time');
        tdim=i;
        
        if strcmp(cfg.toilim, 'all');
            toilim=[1 length(data.time)];
        else
            toilim=[nearest(data.time, cfg.toilim(1)) nearest(data.time, cfg.toilim(end))];
        end
        
        idx_string=[idx_string num2str(toilim(1)) ':' num2str(toilim(end)) ','];
        
    else
        idx_string=[idx_string ':,'];
    end
    
    
    
end

idx_string(end)=[];
dat = data.(cfg.parameter);

dat=eval(['dat(' idx_string ');']);

dat=permute(dat,[chan_dim setxor(chan_dim,1:ndims(dat))]);
datm=nanmean(dat(:,:),2);


% data.(cfg.parameter)=(dat);
% if exist('tdim','var')
%     data.time=data.time(toilim(1):toilim(end));
% end
% if exist('fdim','var')
%     data.freq=data.freq(foilim(1):foilim(2));
% end
% data.dimord='chan';

%%

lay = ft_prepare_layout(cfg, data); 

% import channel positions from layout file
[chNum,X,Y,Width,Height,Lbl] = textread(cfg.layout,'%f %f %f %f %f %s');

% % keep orignal X and Y postions for local topographic plots.  
X0=X;
Y0=Y;
% 
% % normalise channel posiitons to [0 1]
% X=(X-min(X))./(max(X)-min(X));
% Y=(Y-min(Y))./(max(Y)-min(Y));
% 
% % add a border
% b=0.025;
% X=(X-min(X))./(max(X)+(4.*b)-min(X));
% Y=(Y-min(Y))./(max(Y)+(2.*b)-min(Y))+b;
% 
Yd=repmat(Y,1,length(Y))-repmat(Y,1,length(Y))';
Xd=repmat(X,1,length(X))-repmat(X,1,length(X))';
% XYd=sqrt(abs([X Y]*[X Y]' - [[X Y]*[X Y]']'));

XYd=sqrt(sum([Xd(:) Yd(:)].^2,2).*4);
mindist=(min(XYd(find(XYd~=0))));

% X=X./(max(X+mindist+2*b));
% Y=Y./(max(Y+mindist+b));


% do the plotting!

ft_progress('init', 'text')
warning off

datamatrix=NaN.*ones(length(chNum));

for k=1:size(data.labelcmb,1)
 
    
    % identify channelcmb indices for current channl
 %   chancmbidx=chancmbidx+[k.*strcmp(Lbl,data.labelcmb(k,1)) k.*strcmp(Lbl,data.labelcmb(k,2))];
    chanidx1=find(strcmp(Lbl,data.labelcmb(k,1)));
    chanidx2=find(strcmp(Lbl,data.labelcmb(k,2)));
    
    datamatrix(chanidx1,chanidx2)=datm(k);
end

% assume symmetry if matrix is one-sided
if all(isnan(datamatrix(find(triu(ones(size(datamatrix)),1)))));
    datamatrix=tril(datamatrix) + ctranspose(tril(datamatrix));
    if ~isempty(cfg.directionality)
        cfg.directionality=[];
        warning('Connectivity appears to be undirected. Ignoring cfg.directionality');
    end
    
end

% do stuff to deal with complex data
if sum(abs(real(datamatrix(:))))~=0 || sum(abs(imag(datamatrix(:))))~=0
    datamatrix=abs(datamatrix);
elseif sum(abs(real(datamatrix(:))))==0 || sum(abs(imag(datamatrix(:))))~=0
    datamatrix=imag(datamatrix);
end


% get right zlims for all topoplots
if strcmp(cfg.zlim,'maxmin');
    zlim=[min(datm) max(datm)];
elseif strcmp(cfg.zlim,'maxabs');
    zlim=max(abs(datm)).*[-1 1];
elseif strcmp(cfg.zlim,'zeromax');
    zlim=[0 max(datm)];
elseif strcmp(cfg.zlim,'minzero');
    zlim=[min(datm) 0];
end


% go through each channel
for k=1:length(chNum)

if all(isnan(datamatrix(:,k)))
        plotted_chan(k)=0;
else
                plotted_chan(k)=1;
        ft_progress(k/length(chNum), 'Plotting connectivity for sensor %d of %d', k, length(chNum));
        
        %plot the data for current channel
        ft_plot_topo(X0,Y0,datamatrix(:,k),'interpmethod',cfg.interpmethod,...
            'interplim',cfg.interplim,...
            'gridscale',cfg.gridscale,...
            'outline',[],...
            'shading',cfg.shading,...
            'mask',lay.mask,...
            'style',cfg.style,...
            'hpos',X(k),...
            'vpos',Y(k),...
            'width',mindist,...
            'height',mindist);
        


    end
end

ft_progress('close')
warning on
%%
% if colorbar specfied, do one colorbar for all plots
if strcmp(cfg.colorbar,'on');
    colorbar;
end

% change color limits
caxis(cfg.zlim);

%display markers
if ~strcmp(cfg.marker,'off')
  keep_chan = find(plotted_chan);
 
  templay.pos      = [X(keep_chan) Y(keep_chan)];
  templay.width    = mindist;
  templay.height   = mindist;
  templay.label    = Lbl(keep_chan);
  if strcmp(cfg.marker, 'labels') || strcmp(cfg.marker, 'numbers')
    labelflg = 1;
  else
    labelflg = 0;
  end
  if strcmp(cfg.marker, 'numbers')
    for ichan = 1:length(keep_chan)
      templay.label{ichan} = num2str(keep_chan(ichan));
    end
  end
  ft_plot_lay(templay,'box','no','label',labelflg,'point','yes',...
    'pointsymbol',cfg.markersymbol,...
    'pointcolor',cfg.markercolor,...
    'pointsize',cfg.markersize,...
    'labelsize',cfg.markerfontsize,...
    'labeloffset',cfg.labeloffset)
end

axis square;
axis off;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous data
