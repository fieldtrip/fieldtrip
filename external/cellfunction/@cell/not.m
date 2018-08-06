function output = not(input)

output = cellfun(@not, input, 'uniformoutput', false);