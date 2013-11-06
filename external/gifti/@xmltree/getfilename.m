function filename = getfilename(tree)
% XMLTREE/GETFILENAME Get filename method
% FORMAT filename = getfilename(tree)
% 
% tree     - XMLTree object
% filename - XML filename
%__________________________________________________________________________
%
% Return the filename of the XML tree if loaded from disk and an empty 
% string otherwise.
%__________________________________________________________________________
% Copyright (C) 2002-2011  http://www.artefact.tk/

% Guillaume Flandin
% $Id$

filename = tree.filename;
