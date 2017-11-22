function output = sign(input)

output = cellfun(@sign, input, 'uniformoutput', false)