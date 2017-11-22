function output = isreal(input)

output = all(reshape(cell2mat(cellfun(@isreal, input, 'uniformoutput', false)),[],1));
