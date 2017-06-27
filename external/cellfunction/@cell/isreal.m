function output = isreal(input)

output = cellfun(@isreal, input, 'uniformoutput', false);