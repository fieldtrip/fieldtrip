function sens = select_opto(sens, selected)

% SELECT_OPTO make a selection of channels

sens.label = sens.label(selected);
try, sens.chanpos  = sens.chanpos (selected,:); end
try, sens.chanori  = sens.chanori (selected,:); end
try, sens.chantype = sens.chantype(selected,:); end
try, sens.chanunit = sens.chanunit(selected,:); end
try, sens.tra      = sens.tra     (selected,:); end

if ~isfield(sens, 'tra')
  try, sens.optopos  = sens.optopos(selected,:); end
  try, sens.optoori  = sens.optoori(selected,:); end
else
  % tra represents the contribution of each optode to each channel
  % remove the optodes that do not contribute to any channel
  sel      = any(sens.tra~=0,1);
  sens.tra = sens.tra(:,sel);
  try, sens.optopos = sens.optopos(sel,:); end
  try, sens.optoori = sens.optoori(sel,:); end
end

