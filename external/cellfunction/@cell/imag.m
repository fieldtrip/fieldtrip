function output = imag(input)

% output the real part of the data in the input cells
output = cell(size(input));
for k = 1:numel(input)
  output{k} = imag(input{k});
end