function output = times(input1, input2)

if ~iscell(input1)
  input1 = repmat({input1}, size(input2));
end
if ~iscell(input2)
  input2 = repmat({input2}, size(input1));
end
output = cellfun(@times, input1, input2, 'uniformoutput', false);
