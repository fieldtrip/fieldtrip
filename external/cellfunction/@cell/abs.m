function output = abs(input)

output = cellfun(@abs, input, 'uniformoutput', false);