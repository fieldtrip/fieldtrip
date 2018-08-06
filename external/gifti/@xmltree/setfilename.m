function tree = setfilename(tree,filename)
% XMLTREE/SETFILENAME Set filename method
% FORMAT tree = setfilename(tree,filename)
% 
% tree     - XMLTree object
% filename - XML filename
%__________________________________________________________________________
%
% Set the filename linked to the XML tree as filename.
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id: setfilename.m 4460 2011-09-05 14:52:16Z guillaume $

tree.filename = filename;
