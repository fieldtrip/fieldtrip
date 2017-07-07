function output = ctranspose(input)

output = cellfun(@ctranspose, input, 'uniformoutput', false);