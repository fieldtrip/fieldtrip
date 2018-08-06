function neighbours = fixneighbours(neighbours)
% This function converts the old format of the neighbourstructure into the
% new format - although it just works as a wrapper
%
% See also FT_NEIGHBOURSELECTION

neighbours = cellStruct2StructCell(neighbours);

end
