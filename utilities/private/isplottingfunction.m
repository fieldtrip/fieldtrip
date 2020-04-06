function [status] = isplottingfunction(function_name)

% ISPLOTTINGFUNCTION is a helper function for reproducescript, and
% is used for the cfg.reproducescript functionality. It compares the input
% function name with the list of known FieldTrip plotting functions and
% returns 1 if it is a plotting function, and 0 otherwise.

% determine FieldTrip path
ftpath = fileparts(which('ft_defaults'));

functionlist = [];
% browse the main FieldTrip directory
functionlist = [functionlist; dir([ftpath,'/', '*plot*.m'])];
functionlist = [functionlist; dir([ftpath,'/', '*movie*.m'])];
% browse the plotting directory
functionlist = [functionlist; dir([ftpath,'/','plotting/', '*plot*.m'])];
functionlist = [functionlist; dir([ftpath,'/','plotting/', '*movie*.m'])];

% concatenate all function names
functionlist = struct2cell(functionlist);
functionlist(2:end,:) = [];

% compare the function with the list of plotting functions
isinlist = strfind(functionlist, function_name);
status = any(cat(1,isinlist{:}));