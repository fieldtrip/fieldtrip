function sens = undobalancing(sens)

% UNDOBALANCING removes all balancing coefficients from the gradiometer sensor array
%
% This is used in CHANNELPOSITION, FT_PREPARE_LAYOUT

while isfield(sens, 'balance') && isfield(sens.balance, 'current') && ~strcmp(sens.balance.current, 'none')
  fnames = setdiff(fieldnames(sens.balance), 'current');
  indx   = find(ismember(fnames, sens.balance.current));

  if length(indx)==1,
    % undo the synthetic gradient balancing
    fprintf('undoing the %s balancing\n', sens.balance.current);

    % if componentanalysis was followed by rejectcomponent, the balancing matrix is rank deficient
    % leading to problems in the correct allocation of the coils to the channels 
    if strcmp(sens.balance.current, 'invcomp') && strcmp(sens.balance.previous{1}, 'comp')
      tra1 = full(sens.balance.invcomp.tra);
      tra2 = full(sens.balance.comp.tra);
      tmp  = tra1*tra2;
      tmp  = null(tmp); % nullspace after componentanalysis and rejectcomponent
      tmp  = tmp*tmp';  % this is the part which was removed at some point
      [ix,iy] = match_str(sens.balance.comp.labelorg, sens.balance.invcomp.labelnew);
      tra3(iy,iy) = (eye(numel(ix))+tmp(ix,ix))*tra1(iy,iy);
      sens.balance.invcomp.tra = sparse(tra3);
      % FIXME check whether this is robust
    end

    sens = ft_apply_montage(sens, sens.balance.(sens.balance.current), 'inverse', 'yes', 'keepunused', 'yes');
  else
    warning('cannot undo %s balancing\n', sens.balance.current);
    break
  end
end

