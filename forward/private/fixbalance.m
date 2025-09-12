function grad = fixbalance(grad)

% FIXBALANCE ensures that the balancing representation in grad.balance or
% elec.balance field is up to date and consistent, specifically with the
% list of linear projections (or montages) being applied specified as
% cell-array in "current", and not as a string in "current" and cell-array
% in "previous".
%
% See also FT_DATATYPE_SENS

if ~isfield(grad, 'balance')
  grad.balance = struct();
end

if ~isfield(grad.balance, 'current')
  grad.balance.current = {};
end

if ischar(grad.balance.current)
  grad.balance.current = {grad.balance.current};
end

if isfield(grad.balance, 'previous')
  % concatenate them and keep a single list
  grad.balance.current = cat(1, grad.balance.previous(:), grad.balance.current(:))';
  grad.balance = rmfield(grad.balance, 'previous');
end

% remove the 'none'
grad.balance.current = grad.balance.current(~strcmp(grad.balance.current, 'none'));