function mvnx = load_mvnx(filename)
% mvnx = load_mvnx(filename)
% loads a mvnx file
% filename is name of file (including path)
% mvnx is result struct containing all data of the mvnx file
% Xsens Technologies BV 28-05-2015

%% check filename
if isempty(strfind(filename,'mvnx'))
    filename = [filename '.mvnx'];
end
if ~exist(filename,'file')
    error([mfilename ':xsens:filename'],['No file with filename: ' filename ', file is not present or file has wrong format (function only reads .mvnx)'])
end
%% read data
fid = fopen(filename, 'r', 'n', 'UTF-8');
content = char(fread(fid,'char')');
fclose(fid);
%% get data into a cell array, one cell per line
lines = find(content == 10);
cellContent = cell(1,length(lines));

hWaitbar = waitbar(0, 'MVNX import: importing data...');
for n=1:length(lines)
    if n==1
        cellContent{n} = content(1:lines(n)-1);
    else
        cellContent{n} = content(lines(n-1)+1:lines(n)-1);
    end
end

% set file comments on one line
beginComment = find(cellfun(@(x) ~isempty(strfind(x, '<comment>')), cellContent));
endComment = find(cellfun(@(x) ~isempty(strfind(x, '</comment>')), cellContent));

if ~all(beginComment == endComment)
    for i=1:length(beginComment)
        cellContent{beginComment(i)} = [sprintf('%s ',cellContent{beginComment(i):endComment(i)-1}) cellContent{endComment(i)}];
        [cellContent{beginComment(i)+1:endComment(i)}] = deal('<!>');
    end
end
% look for comment lines and remove them
cellContent(cellfun(@(x) x(1)=='-' || (x(1)=='<' && (x(2)=='?' || x(2)=='!')),cellContent)) = [];
% get start of text and clean up some
index = cellfun(@(x) find(x==60,1),cellContent); % lines with only comments on them won't have this mark
openandclose = cellfun(@(x) sum(x==60)==2,cellContent); % lines with have an opening and closing statement in one line
cellContent = cellfun(@(x, y) x(y:end), cellContent, num2cell(index),'UniformOutput',false);

%% get start and end words
n = 0;clear value name
iWord = 0;
word = cell(1,length(cellContent));
wordindex = cell(1,length(cellContent));
wordvalue = cell(1,length(cellContent));
wordfields = cell(1,length(cellContent));
while n < length(cellContent)
    n=n+1;
    line = cellContent{n};oneword = false;
    hooks = find(line == 62);
    iWord = iWord+1;
    wordindex{iWord} = index(n);
    if any(line == 32) && ~openandclose(n)
        if ~isempty(hooks) && hooks(1) < find(line==32,1)
            word{iWord} = line(2:hooks(1)-1);
            iLine = hooks(1)+1;
        else
            word{iWord} = line(2:find(line==32,1)-1);
            iLine = find(line==32,1)+1;
        end
    elseif ~isempty(hooks) && length(line)>=8 && strcmp('comment',line(2:8))
        % add exception for comment
        word{iWord} = line(2:hooks(1)-1);
        iLine = hooks(1)+1;
    elseif openandclose(n)
        word{iWord} = line(2:find(line==62,1)-1);
        iLine = find(line==62,1)+1;
    else
        word{iWord} = line(2:end-1);
        oneword = true;
    end
    if word{iWord}(1) ~= '/'
        if ~oneword && ~openandclose(n)
            k = find(line == 34);
            k = reshape(k,2,length(k)/2)';
            l = [iLine find(line(iLine:end) == 61)+iLine-2];
            fieldname = cell(1,length(l)-1); value = cell(1,length(l)-1);
            if ~isempty(k)
                for il=1:size(k,1)
                    fieldname{il} = line(iLine:find(line(iLine:end) == 61,1)+iLine-2);
                    if size(k,1) > 1 && il < size(k,1)
                        a = strfind(line(iLine:end),'" ')+iLine+1;
                        iLine = a(1);
                    end
                    value{il} = line(k(il,1)+1:k(il,2)-1);
                end
            else
                value = []; fieldname =[];
                value = line(find(line == 62,1)+1:end);
            end
        elseif ~oneword && openandclose(n)
            value = []; fieldname =[];
            value = line(find(line == 62,1)+1:find(line==60,1,'last')-1);
        else
            value = NaN;fieldname = [];
        end
        wordvalue{iWord} = value;
        wordfields{iWord} = fieldname;
    end
end
isendword = cellfun(@(x) x(1) == '/',word);
endwords = cellfun(@(x) x(2:end),word(isendword),'UniformOutput',false);
kindofendwords = unique(endwords);
placeoffirststartword = zeros(1,length(kindofendwords));
placeofallendwords = cell(1,length(kindofendwords));
placeofallstartwords = cell(1,length(kindofendwords));
for n=1:length(kindofendwords)
    lengthWord = length(cell2mat(kindofendwords(n)));
    placeoffirststartword(n) = find(strncmp(word,kindofendwords(n),lengthWord),1);
    placeofallstartwords{n} = find(strncmp(word,kindofendwords(n),lengthWord));
    placeofallendwords{n} = find(strncmp(word,['/' kindofendwords{n}],lengthWord+1));
end
[a b] = sort(placeoffirststartword);
startwords = kindofendwords;
for n=1:length(startwords)
    obj.(startwords{n}).number = length(placeofallstartwords{n});
    obj.(startwords{n}).index = index(a(n));
    obj.(startwords{n}).count = 0;
end

%% get values
for n=1:length(wordvalue)
    if mod(n,5000) == 0
        waitbar((n/length(wordvalue))/2, hWaitbar);
    end
    
    if iscell(wordvalue{n})
        if length(wordvalue{n}) == 1
            B = [];
            try
                B = str2num(wordvalue{n}{1});
            end
            if ~isempty(B)
                wordvalue{n} = B;
            else
                wordvalue{n} = wordvalue{n}{1};
            end
        else
            for m=1:length(wordvalue{n})
                try
                    B = str2num(wordvalue{n}{m});
                    if ~isempty(B)
                        wordvalue{n}{m} = B;
                    end
                end
            end
        end
    else
        try
            B = str2num(wordvalue{n});
            if ~isempty(B)
                wordvalue{n} = B;
            end
        end
    end
end
%% put everything in struct
firstFramesFound = false;
nFramesToBeFound = 5;
nFramesFound = 0;
nFrames = length(cell2mat(placeofallendwords(strcmp(kindofendwords,'frame'))));
lastwordisendword = false;lastendword = '';
superstruct = struct;name = cell(1,max(index)+1);namecounter = ones(1,max(index));endplace = inf(1,max(index));
for iWord=1:length(word)
    if mod(iWord,5000) == 0
        waitbar(iWord/length(word) + 0.5, hWaitbar);
    end
    if word{iWord}(1)=='/'
        startword = word{iWord}(2:end);
        if any(strcmp(startword,name))
            name(find(strcmp(startword,name),1):end) = [];
        end
        if lastwordisendword
            obj.(lastendword).count = 0;
        end
        lastendword = startword;
        lastwordisendword = true;
    else
        lastwordisendword = false;
        try
            if ~isempty(wordfields{iWord})
                wordfields{iWord} = cellfun(@(x) x(x~=32), wordfields{iWord},'UniformOutput',false);
            end
            word(iWord) = checkText(word(iWord));
            wordfields{iWord} = checkText(wordfields{iWord});
            
            if any(strncmp(startwords,word{iWord},length(startwords)))
                if strcmp(word{iWord},'frame')
                    nFramesFound = nFramesFound + 1;
                    if ~firstFramesFound && nFramesFound >= nFramesToBeFound
                        firstFramesFound = true;
                        % preallocate frames in superstruct
                        % Assume remainder of the file contains frames,
                        % with exception of security code at the end
                        frames(1:nFrames,1) = superstruct.mvnx.subject.frames.frame(end);
                        frames(1:nFramesFound-1) = superstruct.mvnx.subject.frames.frame(1:nFramesFound-1);
                    end
                    if firstFramesFound
                        fields = wordfields{iWord};
                        for il = 1:length(fields)
                            frames(nFramesFound).(fields{il}) = wordvalue{iWord}{il};
                        end
                    end
                end
                if isfield(obj, word{iWord})
                    obj.(word{iWord}).count = obj.(word{iWord}).count + 1;
                else
                    obj.(word{iWord}).count = 1;
                end
                name{index(iWord)} = word{iWord};
                namecounter(index(iWord)) = obj.(word{iWord}).count;
                if ~firstFramesFound && (iscell(wordvalue{iWord}) || ~isempty(wordvalue{iWord}) && all(~isnan(wordvalue{iWord})))
                    superstruct = setvalue(superstruct,index(iWord),wordvalue{iWord},name(1:index(iWord)),wordfields{iWord},namecounter);
                end
            elseif iWord > 1 && strcmp(word{iWord},word{iWord-1})
                name{index(iWord)} = word{iWord};
                namecounter(index(iWord)) = namecounter(index(iWord))+1;
                superstruct = setvalue(superstruct,index(iWord),wordvalue{iWord},name(1:index(iWord)),wordfields{iWord},namecounter);
            else
                name{index(iWord)} = word{iWord};
                namecounter(index(iWord)) = 1;
                if firstFramesFound && ~strcmp(word{iWord},'securityCode')
                    frames(nFramesFound).(word{iWord}) = wordvalue{iWord};
                else
                    superstruct = setvalue(superstruct,index(iWord),wordvalue{iWord},name(1:index(iWord)),wordfields{iWord},namecounter);
                end
            end
        catch e
            disp(getReport(e))
        end
    end
end
if firstFramesFound
    superstruct.mvnx.subject.frames.frame = frames';
    if obj.frame.number > nFramesFound
        superstruct.mvnx.subject.frames.frame(nFramesFound+1:end) = [];
    end
end

%% get output
mvnx = superstruct.mvnx;

delete(hWaitbar);
end
function superstruct = setvalue(superstruct,index,value,name,fields,namecounter)
%superstruct = setvalue(superstruct,index,value,name,fields,namecounter)
% add values to a very big struct
if ~isempty(fields) && (~iscell(fields) || ~isempty(fields{1}))
    if length(fields) == 1
        name(end+1) = fields;
    else
        name{end+1} = fields;
    end
    index = index+1;
end
switch index
    case 1
        if ~iscell(name{1})
            superstruct.(name{1}) = value;
        else
            for il = 1:length(name{1})
                superstruct.(name{1}{il}) = value{il};
            end
        end
    case 2
        if ~iscell(name{2})
            superstruct.(name{1})(namecounter(1)).(name{2}) = value;
        else
            for il = 1:length(name{2})
                superstruct.(name{1})(namecounter(1)).(name{2}{il}) = value{il};
            end
        end
    case 3
        if ~iscell(name{3})
            superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})...
                = value;
        else
            for il = 1:length(name{3})
                superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3}{il}) = value{il};
            end
        end
    case 4
        if ~iscell(name{4})
            superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})(namecounter(3))....
                .(name{4}) = value;
        else
            for il = 1:length(name{4})
                superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})(namecounter(3))...
                    .(name{4}{il}) = value{il};
            end
        end
    case 5
        if ~iscell(name{5})
            superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})(namecounter(3))....
                .(name{4})(namecounter(4)).(name{5}) = value;
        else
            for il = 1:length(name{5})
                superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})(namecounter(3))....
                    .(name{4})(namecounter(4)).(name{5}{il}) = value{il};
            end
        end
    case 6
        if ~iscell(name{6})
            superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})(namecounter(3))....
                .(name{4})(namecounter(4)).(name{5})(namecounter(5)).(name{6}) = value;
        else
            for il = 1:length(name{6})
                superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})(namecounter(3))....
                    .(name{4})(namecounter(4)).(name{5})(namecounter(5)).(name{6}{il}) = value{il};
            end
        end
    case 7
        if ~iscell(name{7})
            superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})(namecounter(3))....
                .(name{4})(namecounter(4)).(name{5})(namecounter(5)).(name{6})(namecounter(6))...
                .(name{7}) = value;
        else
            for il = 1:length(name{7})
                superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})(namecounter(3))....
                    .(name{4})(namecounter(4)).(name{5})(namecounter(5)).(name{6})(namecounter(6))...
                    .(name{7}{il}) = value{il};
            end
        end
    case 8
        if ~iscell(name{8})
            superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})(namecounter(3))....
                .(name{4})(namecounter(4)).(name{5})(namecounter(5)).(name{6})(namecounter(6))...
                .(name{7})(namecounter(7)).(name{8}) = value;
        else
            for il = 1:length(name{7})
                superstruct.(name{1})(namecounter(1)).(name{2})(namecounter(2)).(name{3})(namecounter(3))....
                    .(name{4})(namecounter(4)).(name{5})(namecounter(5)).(name{6})(namecounter(6))...
                    .(name{7})(namecounter(7)).(name{8}{il}) = value{il};
            end
        end
end
end

function word = checkText(word)
%word = checkText(word)
% make sure the field name is just text and allowed symbols.
if ~isempty(word)
    if ~iscell(word)
        error('should get cells')
    end
    for i=1:length(word)
        if length(word{i})>1
            word{i}(word{i}=='!') = [];
            word{i}(word{i}==':') = '_';
            if word{i}(end) == '/'
                word{i}(end) = [];
            end
            while (word{i}(1) < 65 || word{i}(1) > 122 || (word{i}(1)> 90 && word{i}(1) < 97))
                word{i}(1)=[];
            end
        end
    end
end
end