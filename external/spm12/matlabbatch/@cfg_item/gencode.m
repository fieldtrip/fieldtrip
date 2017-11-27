function [str, tag, cind] = gencode(item, tag, tagctx)

% function [str, tag, cind] = gencode(item, tag, tagctx)
% Generate code to recreate a cfg_item object. This function calls the
% cfg_item specific gencode_item method. Instead of creating one large,
% deeply nested struct/cell array, it creates separate variables for each
% cfg_item.
%
% Input arguments:
% item - MATLAB variable to generate code for (the variable itself, not its
%        name)
% tag     - optional: name of the variable, i.e. what will be displayed left
%           of the '=' sign. This can also be a valid struct/cell array
%           reference, like 'x(2).y'. If not provided, inputname(1) will be
%           used.
% tagctx  - optional: variable names not to be used (e.g. keywords,
%           reserved variables). A cell array of strings.
%
% Output arguments:
% str  - cellstr containing code lines to reproduce 
% tag  - name of the generated variable
% cind - index into str to the line where the variable assignment is coded
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: gencode.m 3355 2009-09-04 09:37:35Z volkmar $

rev = '$Rev: 3355 $'; %#ok

if nargin < 2
    tag = inputname(1);
end;
if nargin < 3
    tagctx = {};
end
stoptag = tag;
tropts = cfg_tropts({{}},1,inf,1,inf,true);
[str, tag, cind] = gencode_item(item, tag, tagctx, stoptag, tropts);
