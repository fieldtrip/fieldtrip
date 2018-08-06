function res = save(this)
% save an meeg object into a file
% FORMAT res = save(this)
%
% Converts an meeg object to struct and saves it.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: save.m 2394 2008-10-23 15:38:38Z vladimir $


D = struct(this);
D = rmfield(D, 'cache');

res = 1;

try
    save(fullfile(D.path, D.fname), 'D');
catch
    [filename, pathname] = uiputfile('*.mat', 'Select a file to save');
    try
        save(fullfile(pathname, filename), 'D', '-V6');
    catch
        res = 0;
    end
end
