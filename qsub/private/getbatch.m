function [retval] = getbatch

% GETBATCH returns a number that can be used to distinguish the subsequent
% batches. The number automatically increases over subsequent calls.

persistent batch
if isempty(batch)
  batch = 0;
end
batch = batch + 1;
retval = batch;
