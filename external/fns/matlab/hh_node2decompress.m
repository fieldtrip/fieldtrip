function decompress = hh_node2decompress(gridlocs,node_sizes)
% This routine convert the node grid location to the decompress
% vector
% $author: Hung Dang$
% $date: July 21, 2010$
%
% $Log: hh_node2decompress.m$
% Revision 1.1 Wed Jul 21 16:28:47 MDT 2010 hungptit
% New implementation and has been checked.

% The value of the decompress vector is shiffted by 1 to covnert
% from C to MATLAB indexes. 
decompress = gridlocs(:,1) + gridlocs(:,2) * node_sizes(1) + ...
    gridlocs(:,3) * node_sizes(2) * node_sizes(1) + 1;
