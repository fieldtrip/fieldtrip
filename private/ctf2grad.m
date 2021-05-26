function [grad, elec] = ctf2grad(hdr, dewar, coilaccuracy)

% CTF2GRAD converts a CTF header to a gradiometer structure that can be understood by
% the FieldTrip low-level forward and inverse routines. The fieldtrip/fileio
% read_header function can use three different implementations of the low-level code
% for CTF data. Each of these implementations is dealt with here.
%
% Use as
%   [grad, elec] = ctf2grad(hdr, dewar, coilaccuracy)
% where
%   dewar        = boolean, whether to return it in dewar or head coordinates (default is head coordinates)
%   coilaccuracy = empty or a number (default is empty)
%
% See also BTI2GRAD, FIF2GRAD, MNE2GRAD, ITAB2GRAD, YOKOGAWA2GRAD,
% FT_READ_SENS, FT_READ_HEADER

% Copyright (C) 2004-2017, Robert Oostenveld
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

if nargin<2 || isempty(dewar)
  dewar = false;
end

if nargin<3 || isempty(coilaccuracy)
  % if empty it will use the original code
  % otherwise it will use the specified accuracy coil definition from the MNE coil_def.dat
  coilaccuracy = [];
end

if isfield(hdr, 'orig')
  hdr = hdr.orig; % use the original CTF header, not the FieldTrip header
end

% start with empty gradiometer structure
grad = [];
grad.coilpos  = [];
grad.coilori  = [];
grad.tra      = [];

% start with empty electrode structure
elec = [];

% meg channels are 5, refmag 0, refgrad 1, ADCs 18
% UPPT001 is 11
% UTRG001 is 11
% SCLK01  is 17
% STIM    is 11
% SCLK01  is 17
% EEG057  is 9
% ADC06   is 18
% ADC07   is 18
% ADC16   is 18
% V0      is 15

if ~isempty(coilaccuracy)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % use the coil definitions from the MNE coil_def.dat file
  % these allow for varying accuracy which is specified by
  % coilaccuracy = 0, 1 or 2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ft_hastoolbox('mne', 1);
  [ftver, ftpath] = ft_version;
  def = mne_load_coil_def(fullfile(ftpath, 'external', 'mne', 'coil_def.dat'));
  
  k = 1;
  for i=1:length(hdr.res4.senres)
    switch hdr.res4.senres(i).sensorTypeIndex
      case 5 % 5001
        thisdef = def([def.id]==5001 & [def.accuracy]==coilaccuracy);
      case 0 % 5002
        thisdef = def([def.id]==5002 & [def.accuracy]==coilaccuracy);
      case 1 % 5003
        thisdef = def([def.id]==5003 & [def.accuracy]==coilaccuracy);
      otherwise
        % do not add this as sensor to the gradiometer definition
        thisdef = [];
    end % case
    
    if ~isempty(thisdef)
      % the sensors (i.e. the combination of coils that comprises a channel) is
      % defined at [0 0 0] and with the direction [0 0 1]
      pos0 = [0 0 0];
      ori0 = [0 0 1];
      
      if dewar
        % convert from cm to m
        pos2 = hdr.res4.senres(i).pos0(:,1)' / 100;
        if hdr.res4.senres(i).numCoils==2
          % determine the direction from the position of the two coils
          ori2 = (hdr.res4.senres(i).pos0(:,2) - hdr.res4.senres(i).pos0(:,1))';
        else
          % determine the direction from the orientation of the coil
          ori2 = hdr.res4.senres(i).ori0(:,1);
        end
      else
        % convert from cm to m
        pos2 = hdr.res4.senres(i).pos(:,1)' / 100;
        if hdr.res4.senres(i).numCoils==2
          % determine the direction from the position of the two coils
          ori2 = (hdr.res4.senres(i).pos(:,2) - hdr.res4.senres(i).pos(:,1))';
        else
          % determine the direction from the orientation of the coil
          ori2 = hdr.res4.senres(i).ori(:,1);
        end
      end
      ori2 = ori2/norm(ori2);
      
      for j=1:thisdef.num_points
        weight = thisdef.coildefs(j,1);
        pos1 = thisdef.coildefs(j,2:4);
        ori1 = thisdef.coildefs(j,5:7);
        
        [az0,el0,r0] = cart2sph(ori0(1), ori0(2), ori0(3));
        [az2,el2,r2] = cart2sph(ori2(1), ori2(2), ori2(3));
        % determine the homogenous transformation that rotates [1 0 0] towards ori0
        R0 = rotate([0 0 az0*180/pi]) * rotate([0 -el0*180/pi 0]);
        % determine the homogenous transformation that rotates [1 0 0] towards ori2
        R2 = rotate([0 0 az2*180/pi]) * rotate([0 -el2*180/pi 0]);
        % determine the homogenous transformation that rotates ori0 to ori2
        R = R2/R0;
        T = translate(pos2 - pos0);
        H = T*R;
        grad.tra(i,k)     = weight;
        grad.coilpos(k,:) = ft_warp_apply(T*R, pos1); % first the rotation, then the translation
        grad.coilori(k,:) = ft_warp_apply(  R, ori1); % only the rotation
        k = k+1;
      end % for num_points
      
      grad.chanpos(i,:) = pos2;
      grad.chanori(i,:) = ori2;
      % remove the site-specific numbers from each channel name, e.g. 'MZC01-1706' becomes 'MZC01'
      grad.label{i} = strtok(hdr.res4.chanNames(i,:), '-');
    end % if MEG or MEGREF
  end % for each channel
  
  grad.label = grad.label(:);
  grad.unit  = 'm'; % the coil_def.dat file is in meter
  
  remove = cellfun(@isempty, grad.label);
  grad.label   = grad.label(~remove);
  grad.tra     = grad.tra(~remove,:);
  grad.chanpos = grad.chanpos(~remove,:);
  grad.chanori = grad.chanori(~remove,:);
  
  if dewar
    grad.coordsys = 'dewar';
  else
    grad.coordsys = 'ctf';
  end
  
elseif isfield(hdr, 'res4') && isfield(hdr.res4, 'senres')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the header was read using the CTF p-files, i.e. readCTFds
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  sensType  = [hdr.res4.senres.sensorTypeIndex];
  selMEG    = find(sensType==5);
  selREF    = find(sensType==0 | sensType==1);
  selEEG    = find(sensType==9);
  selMEG    = selMEG(:)';
  selREF    = selREF(:)';
  selEEG    = selEEG(:)';
  numMEG    = length(selMEG);
  numREF    = length(selREF);
  numEEG    = length(selEEG);
  
  % determine the number of channels and coils
  coilcount = 0;
  coilcount = coilcount + sum([hdr.res4.senres(selREF).numCoils]);
  coilcount = coilcount + sum([hdr.res4.senres(selMEG).numCoils]);
  chancount = numMEG + numREF;
  % preallocate the memory
  grad.coilpos = zeros(coilcount, 3);         % this will hold the position of each coil
  grad.coilori = zeros(coilcount, 3);         % this will hold the orientation of each coil
  grad.tra     = zeros(chancount, coilcount); % this describes how each coil contributes to each channel
  
  if numEEG>0
    for i=1:numEEG
      n = selEEG(i);
      if dewar
        pos = hdr.res4.senres(n).pos0';
      else
        pos = hdr.res4.senres(n).pos';
      end
      if hdr.res4.senres(n).numCoils~=1
        ft_error('unexpected number of electrode positions in EEG channel');
      end
      % add this position
      elec.elecpos(i       ,:) = pos(1,:);
    end
    % add the electrode names
    elec.label = cellstr(hdr.res4.chanNames(selEEG,:));
  else
    elec = [];
  end
  
  % combine the bottom and top coil of each MEG channel
  for i=1:numMEG
    n = selMEG(i);
    % get coil positions and orientations of this channel (max. 8)
    if dewar
      pos = hdr.res4.senres(n).pos0';
      ori = hdr.res4.senres(n).ori0';
    else
      pos = hdr.res4.senres(n).pos';
      ori = hdr.res4.senres(n).ori';
    end
    if hdr.res4.senres(n).numCoils~=2
      ft_error('unexpected number of coils in MEG channel');
    end
    % add the coils of this channel to the gradiometer array
    grad.coilpos(i       ,:) = pos(1,:);
    grad.coilpos(i+numMEG,:) = pos(2,:);
    grad.coilori(i       ,:) = ori(1,:) .* -sign(hdr.res4.senres(n).properGain);
    grad.coilori(i+numMEG,:) = ori(2,:) .* -sign(hdr.res4.senres(n).properGain);
    grad.tra(i,i       ) = 1;
    grad.tra(i,i+numMEG) = 1;
  end
  
  % the MEG channels always have 2 coils, the reference channels vary in the number of coils
  chancount = 1*numMEG;
  coilcount = 2*numMEG;
  
  % combine the coils of each reference channel
  for i=1:numREF
    n = selREF(i);
    % get coil positions and orientations of this channel (max. 8)
    if dewar
      pos = hdr.res4.senres(n).pos0';
      ori = hdr.res4.senres(n).ori0';
    else
      pos = hdr.res4.senres(n).pos';
      ori = hdr.res4.senres(n).ori';
    end
    % determine the number of coils for this channel
    numcoils = hdr.res4.senres(n).numCoils;
    % add the coils of this channel to the gradiometer array
    chancount = chancount+1;
    for j=1:numcoils
      coilcount = coilcount+1;
      grad.coilpos(coilcount, :)         = pos(j,:);
      grad.coilori(coilcount, :)         = ori(j,:) .* -sign(hdr.res4.senres(n).properGain);
      grad.tra(chancount, coilcount) = 1;
    end
  end
  
  label = cellstr(hdr.res4.chanNames);
  for i=1:numel(label)
    % remove the site-specific numbers from each channel name, e.g. 'MZC01-1706' becomes 'MZC01'
    label{i} = strtok(label{i}, '-');
  end
  
  grad.label = label([selMEG selREF]);
  grad.unit  = 'cm'; % the res4 file represents it in centimeter
  
  if dewar
    grad.coordsys = 'dewar';
  else
    grad.coordsys = 'ctf';
  end
  
  % convert the balancing coefficients into a montage that can be used with the ft_apply_montage function
  if isfield(hdr.BalanceCoefs, 'G1BR')
    meglabel          = label(hdr.BalanceCoefs.G1BR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G1BR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelold  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G1BR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G1BR = montage;
  end
  
  if isfield(hdr.BalanceCoefs, 'G2BR')
    meglabel          = label(hdr.BalanceCoefs.G2BR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G2BR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelold  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G2BR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G2BR = montage;
  end
  
  if isfield(hdr.BalanceCoefs, 'G3BR')
    meglabel          = label(hdr.BalanceCoefs.G3BR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G3BR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelold  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G3BR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G3BR = montage;
  end
  
  if isfield(hdr.BalanceCoefs, 'G3AR')
    meglabel          = label(hdr.BalanceCoefs.G3AR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G3AR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelold  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G3AR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G3AR = montage;
  end
  
  if     all([hdr.res4.senres(selMEG).grad_order_no]==0)
    grad.balance.current = 'none';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==1)
    grad.balance.current = 'G1BR';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==2)
    grad.balance.current = 'G2BR';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==3)
    grad.balance.current = 'G3BR';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==13)
    grad.balance.current = 'G3AR';
  else
    ft_warning('cannot determine balancing of CTF gradiometers');
    grad = rmfield(grad, 'balance');
  end
  
  % sofar the gradiometer definition was the ideal, non-balenced one
  if isfield(grad, 'balance') && ~strcmp(grad.balance.current, 'none')
    % apply the current balancing parameters to the gradiometer definition
    %grad = ft_apply_montage(grad, getfield(grad.balance, grad.balance.current));
    grad = ft_apply_montage(grad, getfield(grad.balance, grad.balance.current), 'keepunused', 'yes');
  end
  
  
elseif isfield(hdr, 'sensType') && isfield(hdr, 'Chan')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the header was read using the open-source MATLAB code that originates from CTF and that was modified by the FCDC
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  selMEG = find(hdr.sensType==5);
  selREF = find(hdr.sensType==0 | hdr.sensType==1);
  selMEG = selMEG(:)';
  selREF = selREF(:)';
  numMEG = length(selMEG);
  numREF = length(selREF);
  
  % combine the bottom and top coil of each MEG channel
  for i=1:numMEG
    n = selMEG(i);
    % get coil positions and orientations of this channel (max. 8)
    if dewar
      pos = cell2mat({hdr.Chan(n).coil.pos}');
      ori = cell2mat({hdr.Chan(n).coil.ori}');
    else
      pos = cell2mat({hdr.Chan(n).coilHC.pos}');
      ori = cell2mat({hdr.Chan(n).coilHC.ori}');
    end
    % determine the number of coils for this channel
    numcoils = sum(sum(pos.^2, 2)~=0);
    if numcoils~=2
      ft_error('unexpected number of coils in MEG channel');
    end
    % add the coils of this channel to the gradiometer array
    grad.coilpos(i       ,:) = pos(1,:);
    grad.coilpos(i+numMEG,:) = pos(2,:);
    grad.coilori(i       ,:) = ori(1,:) .* -sign(hdr.gainV(n));
    grad.coilori(i+numMEG,:) = ori(2,:) .* -sign(hdr.gainV(n));
    grad.tra(i,i       ) = 1;
    grad.tra(i,i+numMEG) = 1;
  end
  numMEGcoils = size(grad.coilpos, 1);
  
  % combine the coils of each reference channel
  for i=1:numREF
    n = selREF(i);
    % get coil positions and orientations of this channel (max. 8)
    if dewar
      pos = cell2mat({hdr.Chan(n).coil.pos}');
      ori = cell2mat({hdr.Chan(n).coil.ori}');
    else
      pos = cell2mat({hdr.Chan(n).coilHC.pos}');
      ori = cell2mat({hdr.Chan(n).coilHC.ori}');
    end
    % determine the number of coils for this channel
    numcoils = sum(sum(pos.^2, 2)~=0);
    % add the coils of this channel to the gradiometer array
    for j=1:numcoils
      grad.coilpos(numMEGcoils+i, :)     = pos(j,:);
      grad.coilori(numMEGcoils+i, :)     = ori(j,:) .* -sign(hdr.gainV(n));
      grad.tra(numMEG+i, numMEGcoils+i) = 1;
    end
  end
  
  grad.label = hdr.label([selMEG selREF]);
  grad.unit  = 'cm'; % the res4 file represents it in centimeter
  
  if dewar
    grad.coordsys = 'dewar';
  else
    grad.coordsys = 'ctf';
  end
  
elseif isfield(hdr, 'sensor') && isfield(hdr.sensor, 'info')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the header was read using the CTF importer from the NIH and Daren Weber
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if dewar
    % this does not work for Daren Webers implementation
    ft_error('cannot return the gradiometer definition in dewar coordinates');
  end
  
  % only work on the MEG channels
  if isfield(hdr.sensor.index, 'meg')
    sel = hdr.sensor.index.meg;
  else
    sel = hdr.sensor.index.meg_sens;
  end
  
  for i=1:length(sel)
    pnt = hdr.sensor.info(sel(i)).location';
    ori = hdr.sensor.info(sel(i)).orientation';
    numcoils(i) = size(pnt,1);
    if size(ori,1)==1 && size(pnt,1)==1
      % one coil position with one orientation: magnetometer
      ori = ori;
    elseif size(ori,1)==1 && size(pnt,1)==2
      % two coil positions with one orientation: first order gradiometer
      % assume that the orientation of the upper coil is opposite to the lower coil
      ori = [ori; -ori];
    else
      ft_error('do not know how to deal with higher order gradiometer hardware')
    end
    
    % add this channels coil positions and orientations
    grad.coilpos = [grad.coilpos; pnt];
    grad.coilori = [grad.coilori; ori];
    grad.label{i} = hdr.sensor.info(sel(i)).label;
    
    % determine the contribution of each coil to each channel's output signal
    if size(pnt,1)==1
      % one coil, assume that the orientation is correct, i.e. the weight is +1
      grad.tra(i,end+1) = 1;
    elseif size(pnt,1)==2
      % two coils, assume that the orientation for each coil is correct,
      % i.e. the weights are +1 and +1
      grad.tra(i,end+1) = 1;
      grad.tra(i,end+1) = 1;
    else
      ft_error('do not know how to deal with higher order gradiometer hardware')
    end
  end
  
  % prefer to have the labels in a column vector
  grad.label = grad.label(:);
  grad.unit  = 'cm'; % the res4 file represents it in centimeter
  
  % reorder the coils, such that the bottom coils are at the first N
  % locations and the top coils at the last N positions. This makes it
  % easier to use a selection of the coils for topographic plotting
  if all(numcoils==2)
    bot = 1:2:sum(numcoils);
    top = 2:2:sum(numcoils);
    grad.coilpos = grad.coilpos([bot top], :);
    grad.coilori = grad.coilori([bot top], :);
    grad.tra = grad.tra(:, [bot top]);
  end
  
else
  ft_error('unknown header to contruct gradiometer definition');
end

% chantype and chanunit are required as of 2016, see FT_DATATYPE_SENS
if ~isempty(grad)
  grad.chantype = ft_chantype(grad);
  grad.chanunit = ft_chanunit(grad);
end
if ~isempty(elec)
  elec.chantype = ft_chantype(elec);
  elec.chanunit = ft_chanunit(elec);
end
