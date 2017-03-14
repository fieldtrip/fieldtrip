function out = floor(in)

out = cellfun(@floor, in, 'uniformoutput', 0);
