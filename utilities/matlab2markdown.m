function matlab2markdown(infile,outfile,varargin)

%This function converts a matlab file to a markdown file, uncommenting all comments,
%assuming they are in correct markdown syntax, plus transforming headings as marked by
%additional % to the markdown # syntax. The script also makes sure code will be highlighted as such.
%As input it takes the name of the file to be converted.
%Best is to provide the full filepath, otherwise it will look for the file within the current path
%The outputfile can be converted back to the original input file using the
%markdown2matlab function.

md_list = '^\*\s|^-\s|^[1-9]*\.\s'; %regex for finding md list syntax
md_indent = '^\s{4}|^\t';

%check input
[inpath, inname, inext]   = fileparts(infile);
if isempty(inpath), [inpath, inname, inext]   = fileparts(which([inname inext])); end
if ~strcmp(inext,'.m')
    error('please specify a matlab file')
end

if nargin < 2
    outfile = infile;
end

%check output
[outpath, outname,outext] = fileparts(outfile);
if isempty(outpath), outpath = inpath; end
if ~strcmp(outext,'.md')
    outext = '.md';
end

%avoid overwriting files
suffix = 1;
newname = outname;
while exist(fullfile(outpath,[newname,outext]),'file') == 2
    newname = [outname sprintf('%0d',suffix)];
    suffix = suffix+1;
end
outname = newname;

infile = fullfile(inpath,[inname,inext]);
outfile = fullfile(outpath,[outname,outext]);

%read & convert file line by line
infid = fopen(infile,'r');
outfid = fopen(outfile,'w');

prevline = '.';
while 1
    tline = fgetl(infid);
    if ~ischar(tline), break, end
    
     %delete included code
        if strcmp(tline,'%include')
            while ~strcmp(tline,'%include_end')
                tline = fgetl(infid);
            end 
            tline = fgetl(infid);
        end
        
        if strcmp(tline,'%output')
            tline = fgetl(infid);
            while ~strcmp(tline,'')
                tline = tline(2:end); %uncomment & indent as if code
                fprintf(outfid,'%s\n',tline);
                tline = fgetl(infid);
            end
        end
    if ~strcmp(tline,'') %ignore blank lines
               
        %each line is either commented out, blank or code
        if any(regexp(tline,md_indent)) && ~any(regexp(strtrim(tline),md_list))
            %ensure that code blocks are lead by a blank line
            if ~any(regexp(prevline,md_indent)) && ~strcmp(prevline,'')
                fprintf(outfid,'%s\n','');
            end
        elseif startsWith(tline,'%')
            %ensure that code blocks are followed by a blank line
            if any(regexp(prevline,md_indent)) && ~any(regexp(strtrim(prevline),md_list))
                fprintf(outfid,'%s\n','');end
            
            tline = tline (2:end);
            %convert headings
            if startsWith(tline,'%')
                indx = strfind(tline,'%');
                while length(indx) > indx(end)
                    indx = indx(1:end-1);
                end
                tline(indx) = strrep(tline(indx),'%','#');
            end
        end
    end
    %%
    fprintf(outfid,'%s\n',tline);
    prevline = tline;
end
fclose(infid);
fclose(outfid);
end