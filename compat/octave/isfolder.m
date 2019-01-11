function tf = isfolder(dirpath)

%ISFOLDER Determine if the input path points to a folder
%   TF = ISFOLDER(PATH) returns true if PATH points to a folder and false otherwise.
%
% This is a compatibility function that should only be added to the path on
% MATLAB versions prior to 2017b.

tf = exist(dirpath,'dir') == 7;

