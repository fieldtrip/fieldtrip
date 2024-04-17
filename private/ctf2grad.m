function [grad, elec] = ctf2grad(hdr, dewar, coilaccuracy, coildeffile)

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
%   coildeffile  = empty or a filename of a valid coil_def.dat file
%
% See also BTI2GRAD, FIF2GRAD, MNE2GRAD, ITAB2GRAD, YOKOGAWA2GRAD,
% FT_READ_SENS, FT_READ_HEADER

% Copyright (C) 2004-2024, Robert Oostenveld
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

if nargin<4 || isempty(coildeffile)
  % if coilaccuracy is not empty, it will use the default coil_def.dat
  % file, which is in external/mne
  coildeffile = [];
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

if isfield(hdr, 'res4') && ~isempty(coilaccuracy)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the header was read using the CTF p-files, i.e. readCTFds
  %
  % use in additionthe coil definitions from the MNE coil_def.dat file
  % these allow for varying accuracy which is specified by coilaccuracy = 0, 1 or 2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ft_hastoolbox('mne', 1);
  if isempty(coildeffile)
    [ftver, ftpath] = ft_version;
    coildeffile     = fullfile(ftpath, 'external', 'mne', 'coil_def.dat');
  end
  def = mne_load_coil_def(coildeffile);

  k = 1;
  for i=1:length(hdr.res4.senres)

    thissens = hdr.res4.senres(i);
    switch thissens.sensorTypeIndex
      case 5 % 5001
        thisdef = def([def.id]==5001 & [def.accuracy]==coilaccuracy);
      case 0 % 5002
        thisdef = def([def.id]==5002 & [def.accuracy]==coilaccuracy);
      case 1 % 5003 or 5004
        % this is a reference gradiometer
        delta = thissens.pos0*[-1;1];

        % determine whether the gradient is 'off' diagonal, comparing the
        % (cosine) of the angle between the line connecting the coils and
        % the orientation of the first coil
        if abs((delta./norm(delta))'*thissens.ori0(:,1))<1e-2
          % it's an off diagonal channel
          thisdef = def([def.id]==5004 & [def.accuracy]==coilaccuracy);
        else
          thisdef = def([def.id]==5003 & [def.accuracy]==coilaccuracy);
        end
      otherwise
        % do not add this as sensor to the gradiometer definition
        thisdef = [];
    end % case

    if ~isempty(thisdef)
      if dewar
        pos = thissens.pos0;
        ori = thissens.ori0;
      else
        pos = thissens.pos;
        ori = thissens.ori;
      end
      pos   = pos./100;   % convert from cm to m

      if thissens.numCoils==2 && thisdef.id~=5004
        % determine the direction from the relative position of the two coils
        %ez = pos(:,2)-pos(:,1);
        %ez = ez./norm(ez);
        ez  = -ori(:,1).*sign(thissens.properGain);

        pos = pos(:,1)'; % take the first coil as local origin
        [ex, ey] = plane_unitvectors(ez);

      elseif thissens.numCoils==2 && thisdef.id==5004
        % 'off' diagonal gradiometer, needs to be treated differently.

        % the local origin is the average of the two coils, and the local
        % x-axis connects the coils
        ex  = pos*[1;-1]; % line between the coils
        pos = mean(pos,2)';

        ez  = -ori(:,1).*sign(thissens.properGain);
        ex  = ex./norm(ex);
        ey  = cross(ez,ex);

        %thisdef.coildefs(:,2:4) = -thisdef.coildefs(:,2:4); % FIXME not sure about this
      else
        % magnetometer coil
        pos = pos';

        % determine the direction from the orientation of the magnetometer coil
        ez = -ori(:,1).*sign(thissens.properGain);
        [ex, ey] = plane_unitvectors(ez);

      end

      for j=1:thisdef.num_points
        weight = thisdef.coildefs(j,1);
        pos1 = thisdef.coildefs(j,2:4);
        ori1 = thisdef.coildefs(j,5:7);

        R = [ex(:) ey(:) ez(:) zeros(3,1);0 0 0 1];
        T = translate(pos);
        grad.tra(i,k)     = weight;
        grad.coilpos(k,:) = ft_warp_apply(T*R, pos1); % first the rotation, then the translation
        grad.coilori(k,:) = ft_warp_apply(  R, ori1); % only the rotation
        k = k+1;
      end % for num_points

      grad.chanpos(i,:) = pos;
      grad.chanori(i,:) = ez;
      grad.label{i}     = hdr.res4.chanNames(i,:);

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
  selMEG    = find(sensType==4 | sensType==5 | sensType==6 | sensType==7);
  selREF    = find(sensType==0 | sensType==1 | sensType==2 | sensType==3);
  selEEG    = find(sensType==9);
  selMEG    = selMEG(:)';
  selREF    = selREF(:)';
  selEEG    = selEEG(:)';
  numMEG    = length(selMEG);
  numREF    = length(selREF);
  numEEG    = length(selEEG);

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

  % determine the number of channels and coils
  chancount = numMEG + numREF;
  coilcount = 0;
  coilcount = coilcount + sum([hdr.res4.senres(selMEG).numCoils]);
  coilcount = coilcount + sum([hdr.res4.senres(selREF).numCoils]);
  % preallocate the memory
  grad.coilpos = zeros(coilcount, 3);         % this will hold the position of each coil
  grad.coilori = zeros(coilcount, 3);         % this will hold the orientation of each coil
  grad.tra     = zeros(chancount, coilcount); % this describes how each coil contributes to each channel

  % keep track of the channels and coils
  chancount = 0;
  coilcount = 0;

  % combine the coils of each MEG channel located in the head shell
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

  grad.label = label([selMEG selREF]);
  grad.unit  = 'cm'; % the res4 file represents it in centimeter
  if ~isempty(elec)
    elec.unit  = 'cm';
  end

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

function [ex, ey] = plane_unitvectors(ez)

% subfunction to obtain a pair of vectors that span the plane orthogonal to
% ez. The heuristic is inspired by the MNE-python code

if abs(abs(ez(3))-1)<1e-5
  ex = [1 0 0]';
else
  ex = zeros(3,1);
  if ez(2)<ez(3)
    if ez(1)<ez(2)
      ex(1) = 1;
    else
      ex(2) = 1;
    end
  else
    if ez(1)<ez(2)
      ex(1) = 1;
    else
      ex(3) = 1;
    end
  end
end
ex = ex - (ex'*ez).*ez;
ex = ex / norm(ex);
ey = cross(ez, ex);
