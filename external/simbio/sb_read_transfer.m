function L = sb_read_transfer(filename)
%reads transfer matrix
mat_filename = [filename,'.mat'];
fid = fopen(mat_filename);
InputText =textscan(fid,'%s','delimiter','\n');
InputText = InputText{1,1};
NumberRows = str2num(strtrim(strrep(InputText{find(strncmp(InputText,'Rows:',5))},'Rows:','')));
NumberColumns = str2num(strtrim(strrep(InputText{find(strncmp(InputText,'Cols:',5))},'Cols:','')));
mab_filename = strtrim(strrep(InputText{find(strncmp(InputText,'FileMatrix:',11))},'FileMatrix:',''));
L.checksum = strtrim(strrep(InputText{find(strncmp(InputText,'FEMchecksum:',12))},'FEMchecksum:',''));
fclose(fid);

fid = fopen(mab_filename);
[L.mat count] = fread(fid, [NumberColumns NumberRows],'double');
L.mat = L.mat';
if(count ~= NumberRows * NumberColumns)
    error('Number of elements read does not match the number of elements expected')
end
fclose(fid);

end
