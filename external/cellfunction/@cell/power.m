function output = power(input, p)

output = cellfun(@power, input, repmat({p}, size(input)), 'uniformoutput', false);