function balance = fixbalance(balance)

% FIXBALANCE ensures that the grad.balance or elec.balance field is
% consistent, specifically with the list of linear projections (or
% montages) being applied specified as cell-array in "current", and not as
% a string in "current" and a cell-array in "previous".
%
% See also FT_DATATYPE_SENS

if ~isfield(balance, 'current')
  balance.current = {};
end

if ischar(balance.current)
  balance.current = {balance.current};
end

if isfield(balance, 'previous')
  % concatenate them and keep a single list
  balance.current = cat(1, balance.previous(:), balance.current(:))';
  balance = rmfield(balance, 'previous');
end

% remove the 'none'
balance.current = balance.current(~strcmp(balance.current, 'none'));