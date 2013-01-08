function s2 = sortfields(s1)

% SORTFIELDS returns the same structure, but with the fields sorted on
% alphabet.

fn = fieldnames(s1);
fn = sort(fn);
for i=1:numel(fn)
  s2.(fn{i}) = s1.(fn{i});
end
