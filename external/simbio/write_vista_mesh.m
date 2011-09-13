function write_vista_mesh()
% Reads a vista format mesh. This is achieved by a correspondently named compiled
% version.
%
% Use as
% write_vista_mesh(filename,nodes,elements,labels(,tensors));
% 
% filename          the name to be saved (with extension .v)
% nodes             n. of nodes*3 field with the position of the nodes
% elements          n. of elements*8 field with the elements
% labels            n. of elements vector with the elements labels
% tensors           (optional, has to be tested) is a n. of elements*6 field
%                   with the tensor conductivity values in order xx-xy-xz-yy-yz-zz.

% Copyright (C) 2011, Johannes Vorwerk, Cristiano Micheli