function markdown2matlab(infile,outfile,varargin)

%This function converts a markdown file to a matlab file, commenting out
%all text, only leaving the code uncommented. As input it takes the name of
%the file to be converted. Best is to provide the full filepath, otherwise
%it will look for the file within the current path
%In addition to the code in the current markdown file it will also add
%any external markdown files that are connected via the include syntax.
%The outputfile can be converted back to the original input file using the
%matlab2markdown function.


md_list = '^\*\s|^-\s|^[1-9]*\.\s'; %regex for finding md list syntax
md_indent = '^\s{4}|^\t';

include = ft_getopt(varargin,'include',0);

%check input
[inpath, inname, inext]   = fileparts(infile);
if isempty(inpath), [inpath, inname, inext]   = fileparts(which([inname inext])); end
if ~strcmp(inext,'.md')
    error('please specify a markdown file')
end

if nargin < 2
    outfile = infile;
end

%check output
[outpath, outname,outext] = fileparts(outfile);
if isempty(outpath), outpath = inpath; end
if ~strcmp(outext,'.m')
    outext = '.m';
end

%avoid overwriting files
if include
    perm = 'a';
else
    perm = 'w';
    suffix = 1;
    newname = outname;
    while exist([newname,outext],'file') == 2
        newname = [outname sprintf('%0d',suffix)];
        suffix = suffix+1;
    end
    outname = newname;
end


infile = fullfile(inpath,[inname,inext]);
outfile = fullfile(outpath,[outname,outext]);
sprintf('writing converted file to %s',outfile)
%read & convert file line by line
infid = fopen(infile,'r');
outfid = fopen(outfile,perm);

while 1 %convert each line in inputfile % write to output file
    tline = fgetl(infid);
    if ~ischar(tline), break, end
    
    %if include shared code
    if ~isempty(strfind(tline,'include')) && ~isempty(strfind(tline,'.md'))
        pathend = strfind(tline,'.md');
        pathbeg = strfind(tline,' /');
        [~,name,ext]= fileparts(tline(pathbeg+1:pathend+2));
        fprintf(outfid,'%s\n',['%' tline]);
        fprintf(outfid,'%s\n','%include');
        fclose(outfid);
        markdown2matlab(which([name,ext]),outfile,'include',1);
        outfid = fopen(outfile,'a');
        fprintf(outfid,'%s\n','%include_end');
        continue
    end
    
    %leave blank lines as is
    if ~strcmp(tline,'')
        %%
        %convert headings
        if startsWith(tline,'#')
            indx = strfind(tline,'#');
            while length(indx) > indx(end)
                indx = indx(1:end-1);
            end
            tline(indx) = strrep(tline(indx),'#','%');
        end
        
        %convert code blocks
        if any(regexp(tline,md_indent)) && ~any(regexp(strtrim(tline),md_list))
            
            %deal with struct output
            if endsWith(tline,'=')
                fprintf(outfid,'%s\n','%output');
                while ~strcmp(tline,'')
                    tline = ['%' tline];
                    fprintf(outfid,'%s\n',tline);
                    tline = fgetl(infid);
                end
            end
            
        else      
            %comment out anything else;
            tline = ['%' tline];
        end
    end
    
    fprintf(outfid,'%s\n',tline);
end

fclose(infid);
fclose(outfid);