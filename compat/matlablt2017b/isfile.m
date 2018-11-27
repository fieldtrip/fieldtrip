function tf = isfile(input)

%ISFILE Determine if the input points to a file
%   TF = ISFILE(INPUT) returns true if INPUT points to a file and false otherwise.
%
% This is a compatibility function that should only be added to the path on
% MATLAB versions prior to 2017b. Note that currently this function only allows
% for a single string in the input, as opposed to the isfile function from MATLAB, which allows for multiple inputs

tf = ismember(exist(input),[2 3 4 5 6]);

