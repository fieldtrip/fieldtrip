function output = sqrt(input)

output = cellfun(@sqrt, input, 'uniformoutput', false);