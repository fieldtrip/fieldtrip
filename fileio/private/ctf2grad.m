function [grad] = ctf2grad(hdr, dewar);

% CTF2GRAD converts a CTF header to a gradiometer structure that can be
% understood by FieldTrip and Robert Oostenveld's low-level forward and
% inverse routines. The fieldtrip/fileio read_header function can use three
% different implementations of the low-level code for CTF data. Each of
% these implementations is dealt with here.
%
% See also READ_HEADER, FIF2GRAD, BTI2GRAD, YOKOGAWA2GRAD

% undocumented option: it will return the gradiometer information in dewar
% coordinates if second argument is present and non-zero

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: ctf2grad.m,v $
% Revision 1.4  2009/08/05 08:36:57  roboos
% preallocate the memory that will hold the coil positions, orientations and weights -> significant speedup
%
% Revision 1.3  2009/03/23 21:16:03  roboos
% don't give error if balancing fails, but only warning and remove the balancing information
%
% Revision 1.2  2009/02/04 13:29:03  roboos
% deal with missing BalanceCoefs in the file using try-catch and isfield (e.g. ArtifactMEG.ds)
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.6  2008/06/26 15:32:19  roboos
% balancing shoudl be subtracted from orig channel, hence added minus sign
%
% Revision 1.5  2008/05/20 12:24:14  roboos
% apply optional balancing to grad.tra
%
% Revision 1.4  2008/05/15 15:11:40  roboos
% add balancing coefficients to the gradiometer definition
%
% Revision 1.3  2008/05/15 13:21:13  roboos
% merged nihm impleemntation into this function
% added new implementation based on CTF p-files
%
% Revision 1.2  2007/03/07 08:57:32  roboos
% use the numeric sensor type for MEG and REF instead of hdr.rowMEG and hdr.rowREF
%
% Revision 1.1  2006/08/31 13:32:11  roboos
% moved from fieldtrip to fileio module
%
% Revision 1.2  2005/06/01 07:59:37  roboos
% added second argument which optionally causes the gradient information to be returend in dewar coordinates
%
% Revision 1.1  2005/05/26 09:55:17  roboos
% renamed the fileio/ctf_grad function to ctf2grad and moved it into fieldtrip/private for consistency with other gradiometer construction functions (nimh2grad and fif2grad)
%
% Revision 1.1  2004/07/02 11:33:47  roboos
% new function that creates a more complete gradiometer definition from the res4 header (compared to read_ctf_res4)
%

% My preferred ordering in the grad structure is:
%   1st 151 coils are bottom coils of MEG channels
%   2nd 151 are the top coils of MEG channels
%   following coils belong to reference channels

if nargin<2 || isempty(dewar)
  dewar = 0;
end

if isfield(hdr, 'orig')
  hdr = hdr.orig; % use the original CTF header, not the FieldTrip header
end

% start with empty gradiometer
grad = [];
grad.pnt = [];
grad.ori = [];
grad.tra = [];

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

if isfield(hdr, 'res4') && isfield(hdr.res4, 'senres')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the header was read using the CTF p-files, i.e. readCTFds
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  sensType  = [hdr.res4.senres.sensorTypeIndex];
  selMEG    = find(sensType==5);
  selREF    = find(sensType==0 | sensType==1);
  selMEG    = selMEG(:)';
  selREF    = selREF(:)';
  numMEG    = length(selMEG);
  numREF    = length(selREF);
  
  % determine the number of channels and coils
  coilcount = 0;
  coilcount = coilcount + sum([hdr.res4.senres(selREF).numCoils]);
  coilcount = coilcount + sum([hdr.res4.senres(selMEG).numCoils]);
  chancount = numMEG + numREF;
  % preallocate the memory
  grad.pnt = zeros(coilcount, 3);         % this will hold the position of each coil
  grad.ori = zeros(coilcount, 3);         % this will hold the orientation of each coil
  grad.tra = zeros(chancount, coilcount); % this describes how each coil contributes to each channel

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
      error('unexpected number of coils in MEG channel');
    end
    % add the coils of this channel to the gradiometer array
    grad.pnt(i       ,:) = pos(1,:);
    grad.pnt(i+numMEG,:) = pos(2,:);
    grad.ori(i       ,:) = ori(1,:) .* -sign(hdr.res4.senres(n).properGain);
    grad.ori(i+numMEG,:) = ori(2,:) .* -sign(hdr.res4.senres(n).properGain);
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
      grad.pnt(coilcount, :)         = pos(j,:);
      grad.ori(coilcount, :)         = ori(j,:) .* -sign(hdr.res4.senres(n).properGain);
      grad.tra(chancount, coilcount) = 1;
    end
  end

  label = cellstr(hdr.res4.chanNames);
  for i=1:numel(label)
    % remove the site-specific numbers from each channel name, e.g. 'MZC01-1706' becomes 'MZC01'
    label{i} = strtok(label{i}, '-');
  end

  grad.label = label([selMEG selREF]);
  grad.unit  = 'cm';

  % convert the balancing coefficients into a montage that can be used with the apply_montage function
  if isfield(hdr.BalanceCoefs, 'G1BR')
    meglabel          = label(hdr.BalanceCoefs.G1BR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G1BR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelorg  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G1BR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G1BR = montage;
  end

  if isfield(hdr.BalanceCoefs, 'G2BR')
    meglabel          = label(hdr.BalanceCoefs.G2BR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G2BR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelorg  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G2BR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G2BR = montage;
  end

  if isfield(hdr.BalanceCoefs, 'G3BR')
    meglabel          = label(hdr.BalanceCoefs.G3BR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G3BR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelorg  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G3BR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G3BR = montage;
  end

  if     all([hdr.res4.senres(selMEG).grad_order_no]==0)
    grad.balance.current = 'none';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==1)
    grad.balance.current = 'G1BR';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==2)
    grad.balance.current = 'G2BR';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==3)
    grad.balance.current = 'G3BR';
  else
    warning('cannot determine balancing of CTF gradiometers');
    grad = rmfield(grad, 'balance');
  end
  
  % sofar the gradiometer definition was the ideal, non-balenced one
  if isfield(grad, 'balance') && ~strcmp(grad.balance.current, 'none')
    % apply the current balancing parameters to the gradiometer definition
    grad = apply_montage(grad, getfield(grad.balance, grad.balance.current));
  end


elseif isfield(hdr, 'sensType') && isfield(hdr, 'Chan')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the header was read using the open-source matlab code that originates from CTF and that was modified by the FCDC
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
      error('unexpected number of coils in MEG channel');
    end
    % add the coils of this channel to the gradiometer array
    grad.pnt(i       ,:) = pos(1,:);
    grad.pnt(i+numMEG,:) = pos(2,:);
    grad.ori(i       ,:) = ori(1,:) .* -sign(hdr.gainV(n));
    grad.ori(i+numMEG,:) = ori(2,:) .* -sign(hdr.gainV(n));
    grad.tra(i,i       ) = 1;
    grad.tra(i,i+numMEG) = 1;
  end

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
      grad.pnt(numMEG+i, :)     = pos(j,:);
      grad.ori(numMEG+i, :)     = ori(j,:) .* -sign(hdr.gainV(n));
      grad.tra(numMEG+i, 2*numMEG+i) = 1;
    end
  end

  grad.label = hdr.label([selMEG selREF]);
  grad.unit  = 'cm';


elseif isfield(hdr, 'sensor') && isfield(hdr.sensor, 'info')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the header was read using the CTF importer from the NIH and Daren Weber
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if dewar
    % this does not work for Daren Webers implementation
    error('cannot return the gradiometer definition in dewar coordinates');
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
      error('do not know how to deal with higher order gradiometer hardware')
    end

    % add this channels coil positions and orientations
    grad.pnt = [grad.pnt; pnt];
    grad.ori = [grad.ori; ori];
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
      error('do not know how to deal with higher order gradiometer hardware')
    end
  end

  % prefer to have the labels in a column vector
  grad.label = grad.label(:);

  % reorder the coils, such that the bottom coils are at the first N
  % locations and the top coils at the last N positions. This makes it
  % easier to use a selection of the coils for topographic plotting
  if all(numcoils==2)
    bot = 1:2:sum(numcoils);
    top = 2:2:sum(numcoils);
    grad.pnt = grad.pnt([bot top], :);
    grad.ori = grad.ori([bot top], :);
    grad.tra = grad.tra(:, [bot top]);
  end

else
  error('unknown header to contruct gradiometer definition');
end
