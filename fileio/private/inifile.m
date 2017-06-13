function readsett = inifile(fileName,operation,keys,style)
%readsett = INIFILE(fileName,operation,keys,style)
%   Creates, reads, or writes data from/to ini (ascii) file.
%
%   - fileName:     ini file name
%   - operation:    can be one of the following:
%                   'new'       (rewrites an existing or creates a new, empty file),
%                   'deletekeys'(deletes keys and their values - if they exist),
%                   'read'      (reads values as strings),
%                   'write'     (writes values given as strings),
%   - keys:         cell array of STRINGS; max 5 columns, min
%                   3 columns. Each row has the same number of columns. The columns are:
%                   'section':      section name string (the root is considered if empty or not given)
%                   'subsection':   subsection name string (the root is considered if empty or not given)
%                   'key':          name of the field to write/read from (given as a string).
%                   'value':        (optional) string-value to write to the ini file in case of 'write' operation
%                   'defaultValue': (optional) value that is returned when the key is not found when reading ('read' operation)
%   - style:        'tabbed' writes sections, subsections and keys in a tabbed style
%                   to get a more readable file. The 'plain' style is the
%                   default style. This only affects the keys that will be written/rewritten.
%
%   - readsett:     read setting in the case of the 'read' operation. If
%                   the keys are not found, the default values are returned
%                   as strings (if given in the 5-th column).
%
%   EXAMPLE:
%   Suppose we want a new ini file, test1.ini with 3 fields.
%   We can write them into the file using:
%
%   inifile('test1.ini','new');
%   writeKeys = {'measurement','person','name','Primoz Cermelj';...
%                'measurement','protocol','id','1';...
%                'application','','description','some...'};
%   inifile('test1.ini','write',writeKeys,'plain');
%
%   Later, you can read them out. Additionally, if any of them won't
%   exist, a default value will be returned (if the 5-th column is given as below).
%   
%   readKeys = {'measurement','person','name','','John Doe';...
%               'measurement','protocol','id','','0';...
%               'application','','description','','none'};
%   readSett = inifile('test1.ini','read',readKeys);
%
%
%   NOTES: When the operation is 'new', only the first 2 parameters are
%   required. If the operation is 'write' and the file is empty or does not exist,
%   a new file is created. When writing and if any of the section or subsection or key does not exist,
%   it creates (adds) a new one.
%   Everything but value is NOT case sensitive. Given keys and values
%   will be trimmed (leading and trailing spaces will be removed).
%   Any duplicates (section, subsection, and keys) are ignored. Empty section and/or
%   subsection can be given as an empty string, '', but NOT as an empty matrix, [].
%
%   This function was tested on the win32 platform only but it should
%   also work on Unix/Linux platforms. Since some short-circuit operators
%   are used, at least Matlab 6.5 should be used.

%
%   FREE SOFTWARE - please refer the source
%   Copyright (c) 2003 by Primoz Cermelj
%   First release on 29.01.2003
%   Primoz Cermelj, Slovenia
%   Contact: primoz.cermelj@email.si
%   Download location: http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=2976&objectType=file
%
%   Version: 1.1.0
%   Last revision: 04.02.2004
%
%   Bug reports, questions, etc. can be sent to the e-mail given above.
%
%   This programme is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation; either version 2
%   of the License, or any later version.
%
%   This programme is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%--------------------------------------------------------------------------

%----------------
% INIFILE history
%----------------
%
% [v.1.1.0] 04.02.2004
% - FIX: 'writetext' option removed (there was a bug previously)
%
% [v.1.01b] 19.12.2003
% - NEW: A new concept - multiple keys can now be read, written, or deleted
%        ALL AT ONCE which makes this function much faster. For example, to
%        write 1000 keys, using previous versions it took 157 seconds on a
%        1.5 GHz machine, with this new version it took only 1.83 seconds.
%        In general, the speed improvement is greater when larger number of
%        read/written keys are considered.
% - NEW: The format of the input parameters has changed. See above.
%
% [v.0.97] 19.11.2003
% - NEW: Additional m-function, strtrim, is no longer needed
%
% [v.0.96] 16.10.2003
% - FIX: Detects empty keys
%
% [v.0.95] 04.07.2003
% - NEW: 'deletekey' option/operation added
% - FIX: A major file refinement to obtain a more compact utility -> additional operations can "easily" be implemented
%
% [v.0.91-0.94]
% - FIX: Some minor refinements
%
% [v.0.90] 29.01.2003
% - NEW: First release of this tool
%
%----------------

global NL_CHAR;

% Checks the input arguments
if nargin < 2
    error('Not enough input arguments');
end
if (strcmpi(operation,'read')) || (strcmpi(operation,'deletekeys'))
    if nargin < 3
        error('Not enough input arguments.');
    end
    if ~exist(fileName)
        error(['File ' fileName ' does not exist.']);
    end
    [m,n] = size(keys);
    if n > 5
        error('Keys argument has too many columns');
    end
    for ii=1:m
        if isempty(keys(ii,3)) || ~ischar(keys{ii,3})
            error('Empty or non-char keys are not allowed.');
        end
    end
elseif (strcmpi(operation,'write')) || (strcmpi(operation,'writetext'))
    if nargin < 3
        error('Not enough input arguments');
    end
    [m,n] = size(keys);
    for ii=1:m
        if isempty(keys(ii,3)) || ~ischar(keys{ii,3})
            error('Empty or non-char keys are not allowed.');
        end
    end
elseif (~strcmpi(operation,'new'))
    error(['Unknown inifile operation: ''' operation '''']);
end
if nargin >= 3
    for ii=1:m
        for jj=1:n
            if ~ischar(keys{ii,jj})
                error('All cells from keys must be given as strings, even the empty ones.');
            end
        end
    end
end
if nargin < 4 || isempty(style)
    style = 'plain';
else
    if ~(strcmpi(style,'plain') || strcmpi(style,'tabbed')) || ~ischar(style)
        error('Unsupported style given or style not given as a string');
    end
end

% Sets new-line character (string)
if ispc
    NL_CHAR = '\r\n';
else
    NL_CHAR = '\n';
end


%----------------------------
% CREATES a new, empty file (rewrites an existing one)
%----------------------------
if strcmpi(operation,'new')
    fh = fopen(fileName,'w');
    if fh == -1
        error(['File: ''' fileName ''' can not be (re)created']);
    end
    fclose(fh);
    return

%----------------------------
% READS key-value pairs out
%----------------------------    
elseif (strcmpi(operation,'read'))
    if n < 5
        defaultValues = cellstrings(m,1);
    else
        defaultValues = keys(:,5);
    end
    readsett = defaultValues;
    keysIn = keys(:,1:3);
    [secsExist,subsecsExist,keysExist,readValues,so,eo] = findkeys(fileName,keysIn);    
    ind = find(keysExist);
    if ~isempty(ind)
        readsett(ind) = readValues(ind);
    end   
    return
    
%----------------------------
% WRITES key-value pairs to an existing or non-existing
% file (file can even be empty)
%----------------------------  
elseif (strcmpi(operation,'write'))
    if m < 1
        error('At least one key is needed when writing keys');
    end
    if ~exist(fileName)
        inifile(fileName,'new');
    end
    writekeys(fileName,keys,style);
    return
    
%----------------------------
% DELETES key-value pairs out
%----------------------------        
elseif (strcmpi(operation,'deletekeys'))
    deletekeys(fileName,keys);
    
    
    
else
    error('Unknown operation for INIFILE.');
end      




%--------------------------------------------------
%%%%%%%%%%%%% SUBFUNCTIONS SECTION %%%%%%%%%%%%%%%%
%--------------------------------------------------


%------------------------------------
function [secsExist,subSecsExist,keysExist,values,startOffsets,endOffsets] = findkeys(fileName,keysIn)
% This function parses ini file for keys as given by keysIn. keysIn is a cell
% array of strings having 3 columns; section, subsection and key in each row.
% section and/or subsection can be empty (root section or root subsection)
% but the key can not be empty. The startOffsets and endOffsets are start and
% end bytes that each key occuppies, respectively. If any of the keys doesn't exist,
% startOffset and endOffset for this key are the same. A special case is
% when the key that doesn't exist also corresponds to a non-existing
% section and non-existing subsection. In such a case, the startOffset and
% endOffset have values of -1.

nKeys = size(keysIn,1);         % number of keys
nKeysLocated = 0;               % number of keys located
secsExist = zeros(nKeys,1);     % if section exists (and is non-empty)
subSecsExist = zeros(nKeys,1);  % if subsection...
keysExist = zeros(nKeys,1);     % if key that we are looking for exists
keysLocated = keysExist;        % if the key's position (existing or non-existing) is LOCATED
values = cellstrings(nKeys,1);  % read values of keys (strings)
startOffsets = -ones(nKeys,1);  % start byte-position of the keys
endOffsets = -ones(nKeys,1);    % end byte-position of the keys

keyInd = find(strcmpi(keysIn(:,1),''));  % key indices having [] section (root section)

line = [];
currSection = '';
currSubSection = '';

fh = fopen(fileName,'r');
if fh == -1
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end

try
    %--- Searching for the keys - their values and start and end locations in bytes
    while 1
       
        pos1 = ftell(fh);
        line = fgetl(fh);
        if line == -1               % end of file, exit
            line = [];
            break
        end
        [status,readValue,readKey] = processiniline(line);        
        if (status == 1)            % (new) section found
            % Keys that were found as belonging to any previous section
            % are now assumed as located (because another
            % section is found here which could even be a repeated one)
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) );  
            if length(keyInd)
                keysLocated(keyInd) = 1;
                nKeysLocated = nKeysLocated + length(keyInd);
            end
            currSection = readValue;
            currSubSection = '';
            % Indices to non-located keys belonging to current section
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) );  
            if ~isempty(keyInd)
                secsExist(keyInd) = 1;
            end
            pos2 = ftell(fh);
            startOffsets(keyInd) = pos2+1;
            endOffsets(keyInd) = pos2+1;
        elseif (status == 2)        % (new) subsection found
            % Keys that were found as belonging to any PREVIOUS section
            % and/or subsection are now assumed as located (because another
            % subsection is found here which could even be a repeated one)
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) & ~keysLocated & strcmpi(keysIn(:,2),currSubSection));
            if length(keyInd)
                keysLocated(keyInd) = 1;
                nKeysLocated = nKeysLocated + length(keyInd);
            end
            currSubSection = readValue;
            % Indices to non-located keys belonging to current section and subsection at the same time
            keyInd = find( ~keysLocated & strcmpi(keysIn(:,1),currSection) & ~keysLocated & strcmpi(keysIn(:,2),currSubSection));
            if ~isempty(keyInd)
                subSecsExist(keyInd) = 1;
            end
            pos2 = ftell(fh);
            startOffsets(keyInd) = pos2+1;
            endOffsets(keyInd) = pos2+1;
        elseif (status == 3)        % key found
            if isempty(keyInd)
                continue            % no keys from 'keys' - from section-subsection par currently in
            end
            currKey = readValue;
            pos2 = ftell(fh);       % the last-byte position of the read key  - the total sum of chars read so far
            for ii=1:length(keyInd)
               if strcmpi( keysIn(keyInd(ii),3),readKey ) && ~keysLocated(keyInd(ii))
                   keysExist(keyInd(ii)) = 1;
                   startOffsets(keyInd(ii)) = pos1+1;
                   endOffsets(keyInd(ii)) = pos2;
                   values{keyInd(ii)} = currKey;
                   keysLocated(keyInd(ii)) = 1;
                   nKeysLocated = nKeysLocated + 1;
               else
                   if ~keysLocated(keyInd(ii))
                       startOffsets(keyInd(ii)) = pos2+1;
                       endOffsets(keyInd(ii)) = pos2+1;
                   end
               end
            end
            if nKeysLocated >= nKeys  % if all the keys are located
                break
            end
        else
                                    % general text found (even empty line(s))
        end       
    %--- End searching
    end    
    fclose(fh);
catch
    fclose(fh);
    error(['Error parsing the file for keys: ' fileName ': ' lasterr]);
end
%------------------------------------




%------------------------------------
function writekeys(fileName,keys,style)
% Writes keys to the section and subsection pair
% If any of the keys doesn't exist, a new key is added to
% the end of the section-subsection pair otherwise the key is updated (changed).
% Keys is a 4-column cell array of strings.

global NL_CHAR;

RETURN = sprintf('\r');
NEWLINE = sprintf('\n');

[m,n] = size(keys);
if n < 4
    error('Keys to be written are given in an invalid format.');
end

% Get keys position first using findkeys
keysIn = keys;    
[secsExist,subSecsExist,keysExist,readValues,so,eo] = findkeys(fileName,keys(:,1:3));

% Read the whole file's contents out
fh = fopen(fileName,'r');
if fh == -1
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
try
    dataout = fscanf(fh,'%c');
catch
    fclose(fh);
    error(lasterr);
end
fclose(fh);

%--- Rewriting the file -> writing the refined contents
fh = fopen(fileName,'w');
if fh == -1
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
try
    tab1 = [];
    if strcmpi(style,'tabbed')
        tab1 = sprintf('\t');
    end
    % Proper sorting of keys is cruical at this point in order to avoid
    % inproper key-writing.
    
    % Find keys with -1 offsets - keys with non-existing section AND
    % subsection - keys that will be added to the end of the file   
    fs = length(dataout);       % file size in bytes
    nAddedKeys = 0;
    ind = find(so==-1);
    if ~isempty(ind)        
        so(ind) = (fs+10);      % make sure these keys will come to the end when sorting
        eo(ind) = (fs+10);
        nAddedKeys = length(ind);
    end
    
    % Sort keys according to start- and end-offsets
    [dummy,ind] = sort(so,1);
    so = so(ind);
    eo = eo(ind);
    keysIn = keysIn(ind,:);
    keysExist = keysExist(ind);
    secsExist = secsExist(ind);
    subSecsExist = subSecsExist(ind);
    readValues = readValues(ind);
    values = keysIn(:,4);   

    % Find keys with equal start offset (so) and additionally sort them
    % (locally). These are non-existing keys, including the ones whose
    % section and subsection will also be added.
    nKeys = size(so,1);
    fullInd = 1:nKeys;
    ii = 1;
    while ii < nKeys
        ind = find(so==so(ii));
        if ~isempty(ind) && length(ind) > 1
            n = length(ind);
            from = ind(1);
            to = ind(end);
            tmpKeys = keysIn( ind,: );
            [tmpKeys,ind2] = sortrows( lower(tmpKeys) );
            fullInd(from:to) = ind(ind2);
            ii = ii + n;
        else
            ii = ii + 1;
        end
    end
    
    % Final (re)sorting
    so = so(fullInd);
    eo = eo(fullInd);
    keysIn = keysIn(fullInd,:);
    keysExist = keysExist(fullInd);
    secsExist = secsExist(fullInd);
    subSecsExist = subSecsExist(fullInd);
    readValues = readValues(fullInd);
    values = keysIn(:,4);    
    
    % Refined data - datain
    datain = [];
    
    for ii=1:nKeys      % go through all the keys, existing and non-existing ones
        if ii==1
            from = 1;   % from byte-offset of original data (dataout)
        else
            from = eo(ii-1);
            if keysExist(ii-1)
                from = from + 1;
            end
        end
        to = min(so(ii)-1,fs);  % to byte-offset of original data (dataout)
        
        if ~isempty(dataout)
            datain = [datain dataout(from:to)];    % the lines before the key
        end
        
        if length(datain) && (~(datain(end)==RETURN || datain(end)==NEWLINE))
            datain = [datain, sprintf(NL_CHAR)];
        end

        tab = [];
        if ~keysExist(ii) 
            if ~secsExist(ii) && ~isempty(keysIn(ii,1))
                if ~isempty(keysIn{ii,1})
                    datain = [datain sprintf(['%s' NL_CHAR],['[' keysIn{ii,1} ']'])];
                end                
                % Key-indices with the same section as this, ii-th key (even empty sections are considered)
                ind = find( strcmpi( keysIn(:,1), keysIn(ii,1)) );
                % This section exists at all keys  corresponding to the same section from know on (even the empty ones)
                secsExist(ind) = 1;
            end
            if ~subSecsExist(ii) && ~isempty(keysIn(ii,2))
                if ~isempty( keysIn{ii,2})
                    if secsExist(ii); tab = tab1;  end;
                    datain = [datain sprintf(['%s' NL_CHAR],[tab '{' keysIn{ii,2} '}'])];
                end
                % Key-indices with the same section AND subsection as this, ii-th key (even empty sections and subsections are considered)
                ind = find( strcmpi( keysIn(:,1), keysIn(ii,1)) & strcmpi( keysIn(:,2), keysIn(ii,2)) );
                % This subsection exists at all keys corresponding to the same section and subsection from know on (even the empty ones)
                subSecsExist(ind) = 1;
            end
        end
        if secsExist(ii) && (~isempty(keysIn{ii,1})); tab = tab1;  end;
        if subSecsExist(ii) && (~isempty(keysIn{ii,2})); tab = [tab tab1];  end;
        datain = [datain sprintf(['%s' NL_CHAR],[tab keysIn{ii,3} ' = ' values{ii}])];
    end
    from = eo(ii);
    if keysExist(ii)
        from = from + 1;
    end
    to = length(dataout);
    if from < to
        datain = [datain dataout(from:to)];
    end
    fprintf(fh,'%c',datain);
catch
    fclose(fh);
    error(['Error writing keys to file: ''' fileName ''' : ' lasterr]);
end
fclose(fh);
%------------------------------------



%------------------------------------
function deletekeys(fileName,keys)
% Deletes keys and their values out; keys must have at least 3 columns:
% section, subsection, and key

[m,n] = size(keys);
if n < 3
    error('Keys to be deleted are given in an invalid format.');
end

% Get keys position first
keysIn = keys;    
[secsExist,subSecsExist,keysExist,readValues,so,eo] = findkeys(fileName,keys(:,1:3));

% Read the whole file's contents out
fh = fopen(fileName,'r');
if fh == -1
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
try
    dataout = fscanf(fh,'%c');
catch
    fclose(fh);
    error(lasterr);
end
fclose(fh);

%--- Rewriting the file -> writing the refined contents
fh = fopen(fileName,'w');
if fh == -1
    error(['File: ''' fileName ''' does not exist or can not be opened.']);
end
try
    ind = find(keysExist);
    nExistingKeys = length(ind);
    datain = dataout;
    
    if nExistingKeys
        % Filtering - retain only the existing keys...
        fs = length(dataout);       % file size in bytes
        so = so(ind);
        eo = eo(ind);
        keysIn = keysIn(ind,:);
        % ...and sorting
        [so,ind] = sort(so);
        eo = eo(ind);
        keysIn = keysIn(ind,:);
        
        % Refined data - datain
        datain = [];
        
        for ii=1:nExistingKeys  % go through all the existing keys
            if ii==1
                from = 1;   % from byte-offset of original data (dataout)
            else
                from = eo(ii-1)+1;
            end
            to = so(ii)-1;  % to byte-offset of original data (dataout)
            
            if ~isempty(dataout)
                datain = [datain dataout(from:to)];    % the lines before the key
            end       
        end
        from = eo(ii)+1;
        to = length(dataout);
        if from < to
            datain = [datain dataout(from:to)];
        end
    end
    
    fprintf(fh,'%c',datain);
catch
    fclose(fh);
    error(['Error deleting keys from file: ''' fileName ''' : ' lasterr]);
end
fclose(fh);
%------------------------------------




%------------------------------------
function [status,value,key] = processiniline(line)
% Processes a line read from the ini file and
% returns the following values:
%   - status:   0   => empty line or unknown string
%               1   => section found
%               2   => subsection found
%               3   => key-value pair found
%   - value:    value-string of a key, section, or subsection
%   - key:      key-string

status = 0;
value = [];
key = [];
line = strim(line);                         % removes any leading and trailing spaces
if isempty(line)                            % empty line
    return
end
if (line(1) == '[') && (line(end) == ']')... % section found
        && (length(line) >= 3)
    value = lower(line(2:end-1));
    status = 1;
elseif (line(1) == '{') &&...                % subsection found
       (line(end) == '}') && (length(line) >= 3)
    value = lower(line(2:end-1));
    status = 2;
else    
    pos = findstr(line,'=');
    if ~isempty(pos)                        % key-value pair found
        status = 3;
        key = lower(line(1:pos-1));
        value = line(pos+1:end);
        key = strim(key);                   % removes any leading and trailing spaces
        value = strim(value);               % removes any leading and trailing spaces
        if isempty(key)                     % empty keys are not allowed
            status = 0;
            key = [];
            value = [];
        end
    end
end
%------------------------------------


%------------------------------------
function outstr = strim(str)
% Removes leading and trailing spaces (spaces, tabs, endlines,...)
% from the str string.
if isnumeric(str);
    outstr = str;
    return
end
ind = find( ~isspace(str) );        % indices of the non-space characters in the str    
if isempty(ind)
    outstr = [];        
else
    outstr = str( ind(1):ind(end) );
end



%------------------------------------
function cs = cellstrings(m,n)
% Creates a m x n cell array of empty strings - ''
cs = cell(m,n);
for ii=1:m
    for jj=1:n
        cs{ii,jj} = '';
    end
end
