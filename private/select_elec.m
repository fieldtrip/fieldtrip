function sens = select_elec(sens, selected)

% SELECT_ELEC make a selection of channels

sens.label = sens.label(selected);
try, sens.chanpos  = sens.chanpos (selected,:); end
try, sens.chanori  = sens.chanori (selected,:); end
try, sens.chantype = sens.chantype(selected,:); end
try, sens.chanunit = sens.chanunit(selected,:); end
try, sens.tra      = sens.tra     (selected,:); end

if ~isfield(sens, 'tra')
  try, sens.elecpos  = sens.elecpos(selected,:); end
  try, sens.elecori  = sens.elecori(selected,:); end
else
  % tra represents the contribution of each electrode to each channel
  % remove the electrodes that do not contribute to any channel
  sel      = any(sens.tra~=0,1);
  sens.tra = sens.tra(:,sel);
  try, sens.elecpos = sens.elecpos(sel,:); end
  try, sens.elecori = sens.elecori(sel,:); end 
end
