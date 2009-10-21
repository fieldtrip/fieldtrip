function [grad] = bti2grad(hdr, balanceflag)

% BTI2GRAD converts a 4D header to a gradiometer structure that can be
% understood by FieldTrip and Robert Oostenveld's low-level forward and
% inverse routines. This function only works for headers that have been
% read using the READ_4D_HDR function.
%
% Use as:
%   [hdr]  = read_4d_hdr(filename);
%   [grad] = bti2grad(hdr);
%
% This function only computes the hardware magnetometer
% definition for the 4D system. This function is based on ctf2grad and
% Gavin Paterson's code, which was adapted from Eugene Kronberg's code
%
% See also FIF2GRAD, CTF2GRAD

% Copyright (C) 2008, Jan-Mathijs Schoffelen 
%
% $Log: bti2grad.m,v $
% Revision 1.5  2009/10/07 09:45:50  jansch
% restructured the handling of balancing; balancing for data in which the
% weight table is of type 1 is still disabled, and balancing is applied for
% data with weight tables of type ~= 1
%
% Revision 1.4  2009/04/02 10:13:15  jansch
% disabled balancing for 148-sensor system (weight table version 1) since
% channel order is not known
%
% Revision 1.3  2009/03/26 10:20:33  jansch
% added balancing based on the weight table used during acquisition. note that
% post acquisition computed weights using 4d software are not incorporated in
% the balancing of the gradiometers
%
% Revision 1.2  2009/01/23 16:15:31  roboos
% removed ; after function declaration
%
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.8  2008/10/20 15:16:16  jansch
% removed the explicit sorting of the channels, this could cause problems
% later on. however, the sensors and references are block-wise sorted still
%
% Revision 1.7  2008/05/15 13:20:36  roboos
% updated documentation
%
% Revision 1.6  2008/05/14 10:20:37  jansch
% included tra-computation when inputting 'm4d' and 'xyz' headers
%
% Revision 1.5  2008/05/14 09:17:04  jansch
% included check for orientation in the case of gradiometers
%
% Revision 1.4  2008/05/14 08:02:40  jansch
% transposed grad.tra (was initially incorrect)
%
% Revision 1.3  2008/05/08 11:10:20  jansch
% implementation in analogy with ctf2grad
%

% for backward compatibility issues FIXME check whether anyone actually uses this code
if isfield(hdr, 'Meg_pos'),
  warning('using outdated code. this does not necessarily lead to correct output');

  grad     = [];
  grad.pnt = hdr.Meg_pos;
  grad.ori = hdr.Meg_dir;
  for i=1:size(grad.pnt,1)
    % grad.label{i} = sprintf('MEG%03d', i);
    grad.label{i} = sprintf('A%d', i); % according to BTi convention
  end
  grad.label = grad.label(:);
  grad.tra = sparse(eye(size(grad.pnt,1)));
  
elseif isfield(hdr, 'config'),
  % hdr has been derived from read_4d_hdr
  
  % it seems as if there is information about the channels at 2 places of
  % the original header. 
  % hdr.channel_data contains info about the actual recorded channels 
  % hdr.config.channel_data contains more important info about ALL 
  %   channels such as position and orientation. 
  % hdr.channel_data maps to hdr.config.channel_data through 
  % hdr.channel_data.chan_no
  
  type    = double([hdr.config.channel_data.type]');
  chan_no = double([hdr.config.channel_data.chan_no]');
  name    = {hdr.config.channel_data.name}'; 
  
  selMEG = chan_no(find(type==1));
  selREF = chan_no(find(type==3));
  selMEG = selMEG(:)';
  selREF = selREF(:)';
  numMEG = length(selMEG);
  numREF = length(selREF);
  
  %sort magnetometers and references
  
  %sortrows does not work here
  %for k = 1:length(selMEG)
  %  n(k) = str2num(name{selMEG(k)}(2:end));
  %end
  
  selALL = [selMEG selREF];
  numALL = length(selALL);
  
  totalcoils = 0;
  numcoils   = zeros(length(type),1);
  for i=1:numALL
    numcoils(selALL(i)) = hdr.config.channel_data(selALL(i)).device_data.total_loops;
  end
  totalcoils = sum(numcoils);
  
  % start with empty gradiometer structure
  grad       = [];
  grad.pnt   = zeros(totalcoils, 3);
  grad.ori   = zeros(totalcoils, 3);
  grad.tra   = zeros(numALL, totalcoils);
  grad.label = cell(numALL,1);
  
  cnt = 0;
  for i=1:numMEG
    n   = selMEG(i);
    pos = cat(2,hdr.config.channel_data(n).device_data.loop_data.position)';
    ori = cat(2,hdr.config.channel_data(n).device_data.loop_data.direction)';
    % determine the number of coils for this channel
    if numcoils(n) ~= size(pos,1)
      error('number of coils does not correspond with number of coil positions');
    end
    % add the coils of this channel to the gradiometer array
    grad.tra(i, cnt+1:cnt+numcoils(n)) = 1;
    % check the orientation of the individual coils in the case of a gradiometer
    % and adjust such that the grad.tra and grad.ori are consistent
    if numcoils(n) > 1,
      c   = ori*ori';
      s   = c./sqrt(diag(c)*diag(c)');
      ori(2:end, :) = ori(2:end, :) .* repmat(-sign(s(2:end,1)), [1 3]);
    end
    for k=1:numcoils(n)
      cnt = cnt+1;
      grad.pnt(cnt,   :) = pos(k,:);
      grad.ori(cnt,   :) = ori(k,:);
    end
    grad.label(i)          = name(n); 
  end
  
  % combine the coils of each reference channel if necessary
  for i=1:numREF
    n   = selREF(i);
    pos = cat(2,hdr.config.channel_data(n).device_data.loop_data.position)';
    ori = cat(2,hdr.config.channel_data(n).device_data.loop_data.direction)';
    % determine the number of coils for this channel
    if numcoils(n) ~= size(pos,1)
      error('number of coils does not correspond with number of coil positions');
    end
    % add the coils of this channel to the gradiometer array
    grad.tra(numMEG+i, cnt+1:cnt+numcoils(n)) = 1; %FIXME check whether ori is OK for gradiometers
    %I think this depends on the orientation of the coils: if they point in opposite directions it's OK
    %check ori by determining the cosine of the angle between the orientations, this should be -1
    % check the orientation of the individual coils in the case of a gradiometer
    % and adjust such that the grad.tra and grad.ori are consistent
    if numcoils(n) > 1,
      c   = ori*ori';
      s   = c./sqrt(diag(c)*diag(c)');
      ori(2:end, :) = ori(2:end, :) .* repmat(-sign(s(2:end,1)), [1 3]);
    end
    for k=1:numcoils(n)
      cnt                = cnt+1;
      grad.pnt(cnt,   :) = pos(k,:);
      grad.ori(cnt,   :) = ori(k,:);
    end
    grad.label(numMEG+i) = {hdr.config.channel_data(n).name}; 
  end
  
  grad.unit  = 'm';
  
  %check whether there is anything to balance
  if ~isa(hdr.user_block_data, 'cell')
    for k = 1:length(hdr.user_block_data)
      tmp{k}=hdr.user_block_data(k);
    end
    hdr.user_block_data = tmp;
  end

  %check whether weights have been applied and stored in header
  for k = 1:length(hdr.user_block_data)
    ubtype{k,1} = hdr.user_block_data{k}.hdr.type;
  end
  ubsel   = strmatch('B_weights_used', ubtype);
  
  if ~isempty(ubsel),
    %balance gradiometers
    weights  = hdr.user_block_data{ubsel};
    if hdr.user_block_data{ubsel}.version==1,
      %the user_block does not contain labels to the channels and references
      %warning('the weight table does not contain contain labels to the channels and references: assuming the channel order as they occur in the header and the refchannel order M.A M.aA G.A');
      label    = {hdr.config.channel_data(:).name}';
      meglabel = channelselection('MEG',    label);
      imeg     = match_str(label, meglabel);
      tabreflabel = {'MxA';'MyA';'MzA';'MxaA';'MyaA';'MzaA';'GxxA';'GyyA';'GyxA';'GzxA';'GzyA'}; %FIXME this is hard coded according to a few tests
      reflabel = channelselection('MEGREF', label);
      [dum, order] = match_str(reflabel, tabreflabel);
      weights.dweights = weights.dweights(imeg,:);
      weights.aweights = weights.aweights(imeg,:);
    else
      meglabel = weights.channames;
      reflabel = weights.drefnames;
      order    = 1:length(reflabel);
    end
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelorg  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -weights.dweights(:,order); zeros(nref, nmeg), eye(nref, nref)];
    balance           = struct(weights.position, montage);
    
    %check the version and the input arguments
    if nargin==1,
      balanceflag = hdr.user_block_data{ubsel}.version~=1;
    end
    
    if balanceflag,
      fprintf('applying digital weights in the gradiometer balancing matrix\n');
      grad.balance      = balance;
      grad.balance.current = weights.position;
      grad              = apply_montage(grad, getfield(grad.balance, grad.balance.current));
    else
      fprintf('not applying digital weights in the gradiometer balancing matrix\n');
    end
  end
  
elseif isfield(hdr, 'grad'),
  %hdr has been derived in a different way and grad is already there, possibly without tra
  grad = hdr.grad;
  if ~isfield(grad, 'tra'), 
    grad.tra = sparse(eye(size(grad.pnt,1)));
  end
end
