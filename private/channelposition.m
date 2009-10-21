function [pnt, lab] = channelposition(sens, varargin)

% CHANNELPOSITION
%
% Use as
%   [pos, label] = channelposition(sens, ...)

% Copyright (C) 2009, Robert Oostenveld & Vladimir Litvak
%
% $Log: channelposition.m,v $
% Revision 1.9  2009/08/10 12:33:57  vlalit
% Adding thresholds to ignore small values in the tra also for Neuromag and planar
%
% Revision 1.8  2009/07/08 07:42:50  jansch
% undone previous adjustment. convention now is that sequential balancing steps
% should be recorded in the grad-structure itself; rather than having
% balance.XXX and balance.YYY, it should be balance.XXX_YYY and
% balance.current = 'XXX_YYY'
%
% Revision 1.7  2009/07/07 08:03:29  jansch
% allowing for sequential unbalancing, convention being that the individual
% balancing steps in balance.current are separated by '_', and the chronological
% ordering of the balancing steps are from right to left
%
% Revision 1.6  2009/06/26 17:39:04  vlalit
% Added the possiblity to handle custom montages applied to MEG sensors (for removing
%  spatial confounds). Hopefully there won't be major side effects.
%
% Revision 1.5  2009/06/03 09:49:03  roboos
% change in whitespace
%
% Revision 1.4  2009/04/03 08:14:27  vlalit
% getting rid of the dependence on statistics toolbox I accidentally introduced by using
%  nanmin.
%
% Revision 1.3  2009/03/30 17:55:17  vlalit
% Changed prepare_layout and headmodelplot to use channelposition. Changed the color
%  of sensor markers in headmodelplot to green for consistency with SPM plots.
%
% Revision 1.2  2009/03/26 16:27:17  roboos
% incorporated the sugegstions by Vladimir
%
% Revision 1.1  2009/03/25 09:04:36  roboos
% new function, will be used as helper function in prepare_layout
%

if isfield(sens, 'balance') && ~strcmp(sens.balance.current, 'none')
  fnames = setdiff(fieldnames(sens.balance), 'current');
  indx   = find(ismember(fnames, sens.balance.current));
  
  if length(indx)==1,
    %  undo the synthetic gradient balancing
    fprintf('undoing the %s balancing\n', sens.balance.current);
    sens = apply_montage(sens, getfield(sens.balance, sens.balance.current), 'inverse', 'yes');
    sens.balance.current = 'none';
  else
    warning('cannot undo %s balancing\n', sens.balance.current);
  end
end

switch senstype(sens)    
  case {'ctf151', 'ctf275' 'bti148', 'bti248'}
    % remove the non-MEG channels altogether
    sel = chantype(sens, 'meg');
    sens.label = sens.label(sel);
    sens.tra   = sens.tra(sel,:);

    % subsequently remove the unused coils
    used = any(abs(sens.tra)<0.5, 1);  % allow a little bit of rounding-off error
    sens.pnt = sens.pnt(used,:);
    sens.ori = sens.ori(used,:);
    sens.tra = sens.tra(:,used);

    % compute distances from the center
    dist = sqrt(sum((sens.pnt - repmat(mean(sens.pnt), size(sens.pnt, 1), 1)).^2, 2));

    % put the corresponding distances instead of non-zero tra entries
    dist = (abs(sens.tra)>0.5).*repmat(dist', size(sens.tra, 1), 1);

    % put nans instead of the zero entries
    dist(~dist) = inf;

    % use the matrix to find coils with minimal distance
    [junk, ind] = min(dist, [], 2);

    lab = sens.label;
    pnt = sens.pnt(ind, :);

  case {'ctf151_planar', 'ctf275_planar', 'bti148_planar', 'bti248_planar'}
    % create a list with planar channel names
    chan = {};
    for i=1:length(sens.label)
      if ~isempty(findstr(sens.label{i}, '_dH')) || ~isempty(findstr(sens.label{i}, '_dV'))
        chan{i} = sens.label{i}(1:(end-3));
      end
    end
    chan = unique(chan);
    % find the matching channel-duplets
    ind = [];
    lab = {};
    for i=1:length(chan)
      ch1 =  [chan{i} '_dH'];
      ch2 =  [chan{i} '_dV'];
      sel = match_str(sens.label, {ch1, ch2});
      if length(sel)==2
        ind = [ind; i];
        lab(i,:) = {ch1, ch2};
        meanpnt1 = mean(sens.pnt(abs(sens.tra(sel(1),:))>0.5, :), 1);
        meanpnt2 = mean(sens.pnt(abs(sens.tra(sel(2),:))>0.5, :), 1);
        pnt(i,:) = mean([meanpnt1; meanpnt2], 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);

  case 'neuromag122'
    % find the matching channel-duplets
    ind = [];
    lab = {};
    for i=1:2:140
      % first try MEG channel labels with a space
      ch1 = sprintf('MEG %03d', i);
      ch2 = sprintf('MEG %03d', i+1);
      sel = match_str(sens.label, {ch1, ch2});
      % then try MEG channel labels without a space
      if (length(sel)~=2)
        ch1 = sprintf('MEG%03d', i);
        ch2 = sprintf('MEG%03d', i+1);
        sel = match_str(sens.label, {ch1, ch2});
      end
      % then try to determine the channel locations
      if (length(sel)==2)
        ind = [ind; i];
        lab(i,:) = {ch1, ch2};
        meanpnt1 = mean(sens.pnt(abs(sens.tra(sel(1),:))>0.5,:), 1);
        meanpnt2 = mean(sens.pnt(abs(sens.tra(sel(2),:))>0.5,:), 1);
        pnt(i,:) = mean([meanpnt1; meanpnt2], 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);

  case 'neuromag306'
    % find the matching channel-triplets
    ind = [];
    lab = {};
    for i=1:300
      % first try MEG channel labels with a space
      ch1 = sprintf('MEG %03d1', i);
      ch2 = sprintf('MEG %03d2', i);
      ch3 = sprintf('MEG %03d3', i);
      sel = match_str(sens.label, {ch1, ch2, ch3});
      % the try MEG channels without a space
      if (length(sel)~=3)
        ch1 = sprintf('MEG%03d1', i);
        ch2 = sprintf('MEG%03d2', i);
        ch3 = sprintf('MEG%03d3', i);
        sel = match_str(sens.label, {ch1, ch2, ch3});
      end
      % then try to determine the channel locations
      if (length(sel)==3)
        ind = [ind; i];
        lab(i,:) = {ch1, ch2, ch3};
        meanpnt1 = mean(sens.pnt(abs(sens.tra(sel(1),:))>0.5,:), 1);
        meanpnt2 = mean(sens.pnt(abs(sens.tra(sel(2),:))>0.5,:), 1);
        meanpnt3 = mean(sens.pnt(abs(sens.tra(sel(3),:))>0.5,:), 1);
        pnt(i,:) = mean([meanpnt1; meanpnt2; meanpnt3], 1);
      end
    end
    lab = lab(ind,:);
    pnt = pnt(ind,:);


  otherwise
    % compute the position for each electrode

    if isfield(sens, 'tra')
      % each channel depends on multiple sensors (electrodes or coils)
      % compute a weighted position for the channel
      [nchan, ncoil] = size(sens.tra);
      pnt = zeros(nchan,3);
      for i=1:nchan
        weight = abs(sens.tra(i,:));
        weight = weight ./ norm(weight);
        pnt(i,:) = weight * sens.pnt;
      end
      lab = sens.label;

    else
      % there is one sensor per channel, which means that the channel position
      % is identical to the sensor position
      pnt = sens.pnt;
      lab = sens.label;
    end

end % switch senstype

n   = size(lab,2);
% this is to fix the planar layouts, which cannot be plotted anyway
if n>1 && size(lab, 1)>1 %this is to prevent confusion when lab happens to be a row array
  pnt = repmat(pnt, n, 1);
end

lab = lab(:);
