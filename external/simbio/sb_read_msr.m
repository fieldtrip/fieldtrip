function b_all = sb_read_msr(filename)
% reads in the output file of the simbio FEM calculation and extracts the
% leadfields

fid = fopen(filename,'r');
ValuesTransposed_flag = 1;

while(ValuesTransposed_flag)
    InputText =textscan(fid,'%s',1,'delimiter','\n');
    if(strfind(InputText{1}{1},'NumberPositions=') == 1)
        N_sens = str2num(deblank(strrep(InputText{1}{1},'NumberPositions=','')));
    end
    
    if(strfind(InputText{1}{1},'NumberTimeSteps=') == 1)
        N_sources = str2num(deblank(strrep(InputText{1}{1},'NumberTimeSteps=','')));
    end
    
    if(strfind(InputText{1}{1},'ValuesTransposed') == 1)
        ValuesTransposed_flag = 0;
    end
end

b_all = zeros(N_sens,N_sources);

for i = 1:N_sources
    InputText = textscan(fid,'%s',1,'delimiter','\n');
    colon_pos = strfind(InputText{1}{1},':');
    a = str2num(InputText{1}{1}(colon_pos+1:end));
    b_all(:,i) = a';
end

fclose(fid);
