function output = exp(input)

output = cellfun(@exp, input, 'uniformoutput', false);