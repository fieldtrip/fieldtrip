function eventout = biosig2fieldtripevent(eventin)

% BIOSIG2FIELDTRIPEVENT converts event information from a biosig hdr into
% fieldtrip events

ft_hastoolbox('BIOSIG', 1);
f = which('sopen.m');
[p,f,e] = fileparts(f);

if contains(p, 'external')
  % assume the biosig to be shipped with the fieldtrip distro, in which
  % case we need to look into fullfile(p,'doc')
else
  % assume the biosig to be from the official distro, in which case we need
  % to look into ../doc, relative to p
  tok = tokenize(p, filesep);
  p   = fullfile(tok{1:end-1});
end

try
  EVT = sopen(fullfile(p, 'doc', 'eventcodes.txt'));sclose(EVT);
catch
  warning('failed to load eventcodes');
  eventout = [];
  return;
end

eventout = struct([]);
if isfield(eventin, 'TYP')
  for index = 1:length( eventin.TYP )
    eType = eventin.TYP(index);
    
    % use file in
    % https://sccn.ucsd.edu/bugzilla/show_bug.cgi?id=1387 to test
    % for boundary events
    if eType < 256 && isfield(eventin,'CodeDesc') && eType < length(eventin.CodeDesc)
      eventout(index).type  = eventin.CodeDesc{eType};
      eventout(index).value = eType;
    elseif isfield(EVT, 'EVENT') && isfield(EVT.EVENT,'CodeIndex') && isfield(EVT.EVENT,'CodeDesc')
      if any(EVT.EVENT.CodeIndex==eType)
        eventout(index).type = EVT.EVENT.CodeDesc{EVT.EVENT.CodeIndex==eType};
      else
        eventout(index).type = 'event';
      end
      eventout(index).value = eType;
      
      % FIXME I leave the below in for reference, I don't know what it
      % means, but we might need it in the future
%       if eType == 32766 || eType == 32767
%         eventout(index).edfplustype = eventout(index).type;
%         eventout(index).type = eeg_boundarytype(event);
%       end
    else
      eventout(index).type  = 'event';
      eventout(index).value = eType;
    end
  end
end
if isfield(eventin, 'POS')
  for index = 1:length( eventin.POS )
    eventout(index).sample = eventin.POS(index);
  end
end
if isfield(eventin, 'DUR')
  if any( [ eventin.DUR ] )
    for index = 1:length( eventin.DUR )
      eventout(index).duration = eventin.DUR(index);
    end
  end
end
if isfield(eventin, 'CHN')
  if any( [ eventin.CHN ] )
    for index = 1:length( eventin.CHN )
      eventout(index).chanindex = eventin.CHN(index);
    end
  end
end

eventout = eventout(:);
