function [outputlabel, outputindex] = inputlabel2outputlabel(cfg, freq)

% INPUTLABEL2OUTPUTLABEL is a subfunction which outputs the cell-arrays
% outputlabel and the corresponding outputindex, and defines how the 
% channels in the original data have to be combined, to provide the 
% wished for combination of the channels, as defined in cfg.combinechan
%
% Configuration-options are:
%   cfg.combinechan = 'planar' combines the horizontal and vertical planar-gradients
%                     'pseudomeg' one gradiometer versus the rest
%   TODO: more flexible way of combining, e.g. by providing a cell-array 

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

if ~isfield(cfg, 'combinechan'), cfg.combinechan = 'no'; end;

if strcmp(cfg.combinechan, 'no')
  % the output labels are similar to the input labels
  outputlabel = {};
  outputindex = {};
  for i=1:length(freq.label)
    outputindex{i} = i;
    outputlabel(i) = freq.label(i);
  end
elseif strcmp(cfg.combinechan, 'planar')
  % find the combination of horizontal and vertical channels that should be combined
  planar      = planarchannelset(freq);
  [sel_dH, H] = match_str(freq.label, planar(:,1));  % indices of the horizontal channels
  [sel_dV, V] = match_str(freq.label, planar(:,2));  % indices of the vertical   channels
  lab_dH      = freq.label(sel_dH);
  lab_dV      = freq.label(sel_dV);

  if length(sel_dH)~=length(sel_dV)
    [int,ih,iv] = intersect(H,V);
    sel_dH = sel_dH(ih);
    sel_dV = sel_dV(iv);
    %error('not all planar channel combinations are complete')
  else
    ih  = 1:length(sel_dH);
    iv  = 1:length(sel_dV);
    int = 1:length(ih); 
  end
  
  % find the other channels that are present in the data
  sel_other = setdiff(1:length(freq.label), [sel_dH(:)' sel_dV(:)']);
  lab_other = freq.label(sel_other);
  
  % define the channel names after combining the planar combinations
  % they should be sorted according to the order of the horizontal planar channels in the data
  [dum, sel_planar] = match_str(freq.label, planar(:,1));
  if exist('ih', 'var'), sel_planar = sel_planar(int); end
  lab_comb          = planar(sel_planar,3);
  [dum, sel_comb]   = match_str(planar(H(ih),3),planar(V(iv),3));
  
  outputlabel = {};
  outputindex = {};

  outputlabel = [outputlabel; lab_other];
  outputindex = [outputindex; mat2cell(sel_other', ones(length(sel_other),1), 1)];

  outputlabel = [outputlabel; lab_comb];
  outputindex = [outputindex; mat2cell([sel_dH sel_dV(sel_comb)], ones(length(sel_dH),1), 2)];
  
elseif strcmp(cfg.combinechan, 'pseudomeg'),
  outputlabel = {};
  outputindex = {};
  for i=1:length(freq.label)
    outputlabel{i} = ['all-',freq.label{i}];
    outputindex{i} = [setdiff(1:length(freq.label),i)];
  end
  for i=1:length(freq.label)
    outputlabel{end+1} = freq.label{i};
    outputindex{end+1} = i;
  end
elseif iscell(cfg.combinechan(1)),
  for i=1:length(cfg.combinechan)
    outputindex{i} = match_str(freq.label, cfg.combinechan{i});
    outputindex{i} = outputindex{i}(:)';
    outputlabel{i} = cell2mat(freq.label(outputindex{i})');
  end
else
  error('unknown combination method');
end
