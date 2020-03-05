function display(this)
% Display method for GIfTI objects
% FORMAT display(this)
% this   -  GIfTI object
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: display.m 4182 2011-02-01 12:29:09Z guillaume $

display_name = inputname(1);
if isempty(display_name)
    display_name = 'ans';
end

if length(this) == 1 && ~isempty(this.data)
    eval([display_name ' = struct(this);']);
    eval(['display(' display_name ');']);
else
    disp(' ')
    disp([display_name ' =']);
    disp(' ');
    eval([display_name ' = this;']);
    eval(['disp(' display_name ');']);
end