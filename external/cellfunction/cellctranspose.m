function [y] = cellctranspose(x)

y = cellfun(@ctranspose, x, 'uniformoutput', false);
