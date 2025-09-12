function sens = undobalancing(sens)

% UNDOBALANCING removes all balancing coefficients from the gradiometer sensor array
%
% This is used in CHANNELPOSITION, FT_PREPARE_LAYOUT, FT_SENSTYPE

% ensure that the representation of the balancing is up to date
sens = fixbalance(sens);
fprintf('.')

while ~isempty(sens.balance.current)
  % undo the synthetic gradient balancing
  fprintf('undoing the %s balancing for the gradiometer definition\n', sens.balance.current{end});

  % if ft_componentanalysis was followed by ft_rejectcomponent, the balancing matrix is rank deficient
  % leading to problems in the correct allocation of the coils to the channels
  if length(sens.balance.current)>1 && strcmp(sens.balance.current{end}, 'invcomp') && strcmp(sens.balance.current{end-1}, 'comp')
    tra1 = full(sens.balance.invcomp.tra);
    tra2 = full(sens.balance.comp.tra);
    tra3 = tra1;
    tmp  = tra1*tra2;
    tmp  = null(tmp); % nullspace after ft_componentanalysis and ft_rejectcomponent
    tmp  = tmp*tmp';  % this is the part which was removed at some point
    [ix,iy]     = match_str(sens.balance.comp.labelold, sens.balance.invcomp.labelnew);
    tra3(iy,iy) = (eye(numel(ix))+tmp(ix,ix))*tra1(iy,iy);
    sens.balance.invcomp.tra = tra3;
    % FIXME check whether this is robust
  end

  if strcmp(sens.balance.current{end}, 'planar')
    if isfield(sens, 'type') && ~isempty(strfind(sens.type, '_planar'))
      % remove the _planar postfix from the sensor type
      sens.type = sens.type(1:(end-7));
    end
  end

  sens = ft_apply_montage(sens, ft_inverse_montage(sens.balance.(sens.balance.current{end})), 'keepunused', 'yes', 'warning', 'no');
  sens.balance.current = sens.balance.current(1:end-1);

  if ~isfield(sens, 'chanpos') || any(isnan(sens.chanpos(:)))
    % this happens if the data has been component-analyzed
    % try to reconstruct the channel position and orientation
    [pos, ori, lab] = channelposition(sens);
    [sel1, sel2] = match_str(sens.label, lab);
    sens.chanpos(sel1,:) = pos(sel2,:);
    if isfield(sens, 'chanori')
      % this only applies to meg data
      sens.chanori(sel1,:) = ori(sel2,:);
    end
  end

end % while
