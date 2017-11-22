function output = log(input)

output = cellfun(@log, input, 'uniformoutput', false);