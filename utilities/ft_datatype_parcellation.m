function parcellation = ft_datatype_parcellation(parcellation, varagin)

% FT_DATATYPE_PARCELLATION describes the FieldTrip MATLAB structure for
% parcellated cortical surface-based data and atlases.
%
% A parcellation describes the tissue types for each of the surface
% elements. Parcellations are often, but not always labeled. A parcellatoin
% can be used to estimate the activity from MEG data in a known region of
% interest. A surface-based atlas is basically a very detailled parcellation
% with anatomical labels for each tissue type.
%
% An example of a surface based Brodmann parcellation looks like this
%
%              pos: [8192x3]         vertices of the cortical sheet
%              tri: [16382x3]        triangles of te cortical sheet
%            coord: 'ctf'            the (head) coordinate system in which the vertex positions are expressed
%             unit: 'mm'             the units in which the coordinate system is expressed
%         brodmann: [8192x1 uint8]   values from 1 to N, the value 0 means unknown
%    brodmannlabel: {Nx1 cell}
%
% The only difference to the source data structure is that the parcellation
% structure contains the additional fields XXX. See FT_DATATYPE_SOURCE for
% further details.
