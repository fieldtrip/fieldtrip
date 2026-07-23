function sens = select_grad(sens, selected)

% SELECT_GRAD make a selection of channels

sens.label = sens.label(selected);
try, sens.chanpos  = sens.chanpos (selected,:); end
try, sens.chanori  = sens.chanori (selected,:); end
try, sens.chantype = sens.chantype(selected,:); end
try, sens.chanunit = sens.chanunit(selected,:); end
try, sens.tra      = sens.tra     (selected,:); end

if ~isfield(sens, 'tra')
  try, sens.coilpos  = sens.coilpos(selected,:); end
  try, sens.coilori  = sens.coilori(selected,:); end
else
  % tra represents the contribution of each coil to each channel
  % remove the coils that do not contribute to any channel
  sel      = any(sens.tra~=0,1);
  sens.tra = sens.tra(:,sel);
  try, sens.coilpos = sens.coilpos(sel,:); end
  try, sens.coilori = sens.coilori(sel,:); end
end
