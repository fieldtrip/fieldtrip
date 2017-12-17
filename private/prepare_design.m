function [cfg] = prepare_design(cfg)

% PREPARE_DESIGN makes a design matrix on the basis of the information in
% cfg (i.c., cfg.statistic, cfg.ext, and an initial design in cfg.design) 
% and puts this design matrix in cfg.design. PREPARE_DESIGN also gives default
% values for cfg.ivar, which specifies the independent variable, and cfg.uvar, 
% which specifies the units-of-observation.
%
% PREPARE_DESIGN will be called from STATISTICS_WRAPPER whenever the user
% has not specified the cfg.design field.
%
% To construct the design matrix, PREPARE_DESIGN has to know whether 
% cfg.statistic is a statistic for a between- or a within-units 
% design. This is because, for the calculation of a statistic for a 
% within-units design, the unit-of-observation to which a particular
% replication belongs has to be known. PREPARE_DESIGN determines the design
% type (between or within) on the basis of cfg.statistic.
%
% The design type has implications for how the data have to be passed to 
% PREPARE_DESIGN:
% 1. For a between-units design, by default, cfg.design is equal to the
%   last column of cfg.design. (If cfg.design is produced by
%   PREPARE_TIMEFREQDATA, and the varargin-argument of PREPARE_TIMEFREQDATA 
%   contains one data set for every condition, then this column 
%   contains the rank orders of these data sets.) This default 
%   option can be overruled by cfg.ext, which contains an
%   external variable of order Nreplications X Nextvar. (Nextvar is the number of 
%   external variables, and this can be larger that 1.) The order of the
%   replications is determined by the order of the data sets in varargin: the
%   replications in varargin{1} come first, followed by the replications
%   in varargin{2}, etc. In the case of multiple external variables, by specifying 
%   cfg.ivar and cfg.cvar, the independent and the control variables can be specified.
%
% 2. For a within-units design, the default option is the following: (1)
%   the independent variable is equal to the last column of data.design,
%   and (2) the unit-variable is equal to the next-to-last column of
%   data.design. This default option only makes sense if the
%   varargin-argument of PREPARE_TIMEFREQDATA contains one data set for
%   every condition, and if the units in these data sets (subjects or
%   trials) correspond to each other. This default option can be overruled by 
%   cfg.ext, which has order Nwcond X Nextvar or (Nunits*Nwcond) X Nextvar.
%   (Nwcond is the number of within-unit conditions.) The default option of 
%   comparing all within-units conditions with each other can be overruled by 
%   specifying 'ivar' and 'cvar'.

% Copyright (C) 2006, Eric Maris
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

% determine whether a beween or a within-units design is requested.
if any(strcmp(cfg.statistic,{'indepsamplesT','indepsamplesregrT','indepsamplesZcoh','indepsamplesF'}))
    designtype = 'between';
end
if any(strcmp(cfg.statistic,{'depsamplesregrT','depsamplesT','actvsblT','depsamplesFmultivariate'}))
    designtype = 'within';
end
if ~(exist('designtype')==1)
    ft_warning('Unknown test statistic.');
    return
end

if ~isfield(cfg,'design')
  ft_error('You should specify an initial design in cfg.design.');
end
initialdesign=cfg.design;

if strcmp(designtype,'between')   % between-units conditions
  % construct the design matrix
  cfg.design=initialdesign(:,end); % the default option
  if ~isfield(cfg,'ivar')
    cfg.ivar = 1;
  end
  % overwrite the default option if there is a .ext-field in the
  % configuration
  if isfield(cfg,'ext')
    if size(cfg.ext,1)~=size(initialdesign,1)
        ft_error('Incompatible number of replications in cfg.ext.');
    end
    design=cfg.ext;
  end

elseif strcmp(designtype,'within')  % within-units conditions
  % construct the design matrix
  % start with the default option
  nrepl=size(data.biol,1);
  cfg.design=zeros(nrepl,2);
  cfg.design(:,1)=initialdesign(:,end);
  cfg.design(:,2)=initialdesign(:,(end-1));
  % perform some checks on the design
  wcondlabels=sort(cell2mat(unique(cfg.design(:,1))));
  nwcond=length(wcondlabels);
  unitlabels=sort(cell2mat(unique(cfg.design(:,2))));
  nunits=length(unitlabels);
  for wcondindx=1:nwcond
    selvec = (cfg.design(:,1)==wcondlabels(wcondindx));
    unitsthiscond = sort(cfg.design(selvec,2));
    if length(unitthiscond)==length(unitlabels) && ~all(unitthiscond==unitlabels)
      ft_error('The last two columns of initialdesign do not specify a within-units design.');
    end
  end
  if ~isfield(cfg,'ivar')
    cfg.ivar = 1;
  end
  if ~isfield(cfg,'uvar')
    cfg.uvar = 2;
  end
  % overwrite the default option if there is a .ext-field in the
  % configuration
  if isfield(cfg,'ext')
    dimext=size(cfg.ext);
    cfg.design=zeros(nrepl,dimext(2)+1);
    if dimext(1)==nwcond
      for wcondindx=1:nwcond
        for unitindx=1:nunits
          cfg.design(((wcondindx-1)*nunits + unitindx),[1:dimext(2)])=cfg.ext;
        end
      end
    elseif dimext(1)==nrepl
      cfg.design(:,[1:dimext(2)])=cfg.ext;
    else
      ft_error('The number of rows in cfg.ext must be equal to the number of conditions or the number of replications (number of conditions times number of units-of-observation).');
    end
    cfg.design(:,dimext(2)+1)=initialdesign(:,(end-1));
    cfg.uvar = size(design,2);
  end
end

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

