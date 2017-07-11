function structCell = cellStruct2StructCell(cellStruct)
% Converts a cell array of structure arrays into a structure array

structCell = struct;
for i=1:numel(cellStruct)
    if i==1, structCell = cellStruct{i};    end
    structCell(i) = cellStruct{i}; 
end
