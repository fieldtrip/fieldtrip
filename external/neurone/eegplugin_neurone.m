% eegplugin_neurone() - EEGLAB plugin to import data from a NeurOne device.
% 
% Usage:
%   >> eegplugin_neurone(fig,try_strings,catch_strings)
%
% Inputs:
%   fig            - [integer]  handle to EEGLAB figure
%   try_strings    - [struct] "try" strings for menu callbacks.
%   catch_strings  - [struct] "catch" strings for menu callbacks.
%
% NeurOne data import plugin consists of the following files:
%   pop_readneurone.m
%   guireadneurone.m
%   guireadneurone.fig
%   readneurone.m
%   readneuronedata.m
%   readneuroneevents.m
%   neurone_logo.png
%   mega_gradient_edit.png
%
% This plugin was created according to the instructions provided by the
% creators of EEGLAB. These instructions can be found e.g. from the website:
% http://sccn.ucsd.edu/wiki/A07:_Contributing_to_EEGLAB
% 
% Current version: 1.0.3.4 (2016-06-17)
% Author: Mega Electronics

function vers = eegplugin_neurone(fig,try_strings,catch_strings)

vers='NeurOne data import 1.0.3.4';

% Check the number of input arguments
if nargin < 3
    error('Not enough input arguments.');
end

% Add plugin folder to path
if  ~exist('pop_readneurone.m')               
    path=which('eegplugin_neurone.m');
    [path filename]=fileparts(path);
    path=[path filesep];
    addpath([path version] );
end

% Find the 'Import data' -menu
importmenu=findobj(fig,'tag','import data');

% Construct command
cmd = [try_strings.no_check '[EEG LASTCOM]=pop_readneurone;' catch_strings.new_and_hist];

% Create the menu for NeurOne import
uimenu(importmenu,'label','From a NeurOne file (.ses)',...
    'Callback',cmd,'separator','on');

