function output = isnan(input)

output = cellfun(@isnan, input, 'uniformoutput', false);