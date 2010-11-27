% FastICA for Matlab 7.x and 6.x
% Version 2.5, October 19 2005
% Copyright (c) Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen.
%
% Type fasticag to launch the graphical user interface
%
% Please refer to your Matlab documentation on how to add FastICA to your
% Matlab search path. (One place to start is the path-command)
%
% FastICA programs:
%   FASTICAG  - Graphical user interface for FastICA
%   FASTICA   - command line version of FastICA
%
% Separate functions used by FastICA programs.
%   FPICA     - main algorithm for calculating ICA
%   WHITENV   - function for whitening data
%   PCAMAT    - calculates the PCA for data
%   REMMEAN   - function for removing mean
%
%   GUI_CB    - needed by fasticag
%   GUI_ADV   - needed by fasticag
%   GUI_ADVC  - needed by fasticag
%   GUI_L     - needed by fasticag
%   GUI_LC    - needed by fasticag
%   GUI_S     - needed by fasticag
%   GUI_SC    - needed by fasticag
%   GUI_CG    - needed by fasticag
%   GUI_HELP  - needed by fasticag
%
%   ICAPLOT   - for plotting the signals
%               (also used by fastica and fasticag)
%
% Misc.
%   DEMOSIG   - generates some test signals
%
% Deprecated
%   dispsig   - plots the data vectors
%               replaced by icaplot

% @(#)$Id$