function L = hh_leadfield_new(recipdata,compress,node_sizes,nnode)
% hh_leadfield is a function to extract the lead field matrix from the
% reciprocity data matrix
% Use as
% L = hh_leadfield_new(recipdata,compress,node_sizes,nnode)
% 
% Descriptions:
% 1. L:            Lead field matrix
% 2. recipdata:    Reciprocity data matrix
% 3. compress:     Compress vector
% 4. nodesizes:    The sizes of the node volume
% 5. nnode:        node position in the compress vector
%
% $Id$
% $Date: July 23, 2008$
% Copyright 2008 by Hung Dang        
% $Log: hh_leadfield_new.m$
% Revision 1.1 Wed Jul 21 16:38:48 MDT 2010
% Copy from the old HHSIM routine and update the document. It has been checked with hh_leadfield.m.
    
% Init parameters
ROW = node_sizes(1);
COL = node_sizes(2);

% Update position for the head and tail of three components
head_x = compress(nnode + 1);
tail_x = compress(nnode - 1);
head_y = compress(nnode + ROW);
tail_y = compress(nnode - ROW);
head_z = compress(nnode + COL * ROW);
tail_z = compress(nnode - COL * ROW);

% Extract the lead field matrix
if ((head_x > 0) && (tail_x > 0) && (head_y > 0) && (tail_y > 0) && (head_z > 0) && (tail_z > 0))
    % Extract the lead field matrix for x, y and z component
    L = recipdata(:,[head_x,head_y,head_z]) - recipdata(:,[tail_x,tail_y,tail_z]);
else
    L = [];
end
