function output = isfinite(input)

output = cellfun(@isfinite, input, 'uniformoutput', false);