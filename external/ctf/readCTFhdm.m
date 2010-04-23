function hdm=readCTFhdm(fileName);

%  Version 1.1   19 April 2007 - Test date.
%                21 March 2007.  Modified to read v6.0 .hdm files that have additional 
%                                fields in MultiSphere_Data
%  Reads a head model file and returns the contents as a structure.  The purpose is
%  to make localSpheres head models available in the MATLAB environment.
%  Head Model File format is a "Config Reader" file.  It is defined in document 
%  "CTF MEG FIle Formats', PN900-0088.

%  Input : fileName: Name of a .hdm file

%  Output: hdm: A structure with a field for comments, and a field for each 
%               class in the .hdm file.  Within each class field, there is a field for 
%               each tag found in the class, except for the MultiSphere_Data where 
%               the information is repackaged as SQUIDname, sphereOrigin and radius.

%  NOTE : In ths Config reader file format, white space can be spaces or tabs or 
%  multiple spaces or multiple tabs.  Therefore cannot know where the '{' and '}'
%  characters appear.

persistent printWarning

%  Create a list of characters that are allowed in classnames.
%  his wil be used to decide if the file is not in "Config Reader" format.
allowedChars=char(['A':'Z' 'a':'z' '0':'9' '_']);
forbiddenChars=setdiff(char(0:255),allowedChars);

%  Check inputs and print messages.
if nargin==0 & nargout==0
    fprintf('readCTFhdm   Version 1.1   19 April 2007\n');
    fprintf(['\thdm=readCTFhdm(fileName) creates a structure containing the contents of ',...
            'a CTF-format head model file.\n',...
            '\t\t\tfileName is the name of a head model file including path and ',...
            'extension .hdm\n',...
            '\t\t\thdm is a MATLAB structure containing the head model (comments,source \n',...
            '\t\t\t\tinformation, single-sphere model, multiple local Sphere model).\n',...
            '\t\t\tThe head model file format is described in document ',...
            '"CTF MEG File Formats", PN900-0088.\n\n']);
    return
end

hdm=struct([]);  % Force error in calling program with an early return

if nargout~=1 | nargin~=1
    fprintf(['\nreadCTFhdm:  %d input and %d output arguments.\n',...
            '             readCTFhdm must have 1 input and 1 output.\n\n'],nargin,nargout);
    if nargout==0;clear hdm;end
    return
elseif ~ischar(fileName) | isempty(fileName) | ndims(fileName)>2 | min(size(fileName))~=1
    fprintf('readCTFhdm: Input argument is not character string, or is empty.\n');
    return
elseif exist(fileName)~=2
    fprintf('readCTFhdm: Cannot find file %s\n',fileName);
    return
end

fid=fopen(fileName,'r','ieee-be');

hdm=struct('note',char([]));
classOpen=0;
className=char([]);
while 1
    txtline=fgetl(fid);
    if isequal(txtline,-1);break;end  % eof of file
    txtline=deblank(txtline);
    if isempty(deblank(txtline))    % Empty line
        continue
    end
    % Remove leading white space
    txtline=deblank(txtline(length(txtline):-1:1));
    txtline=txtline(length(txtline):-1:1);
    
    %  Is this line a comment?  If yes it will be added to a note field.
    comment=strcmp('//',sscanf(txtline,'%s',1));
    if comment & classOpen==0;   % Comment
        hdm.note=strvcat(hdm.note,txtline);
    elseif strcmp(txtline,'{')   % Start of new class
        if isempty(className);
            fprintf('readCTFhdm: Error decoding file %s.\n',filename);
            fprintf('\tEncountered ''{'' without first finding a class name.\n');
            hdm=struct([]);  % Force an error in the calling program.
            break
        end
        classOpen=1;
    elseif strcmp(txtline,'}')  % Close the class
        className=char([]);
        classOpen=0;
    elseif classOpen==0   % Create a new class.
        className=txtline;
        if ~strcmp(className,strtok(className,forbiddenChars))
            fprintf('\nreadCTFhdm : File %s\n',fileName);
            fprintf(['             Encountered className ''%s'' which includes characters ',...
                    'that are not allowed in a class name.\n\n'],className);
            break
        end
        % Create an empty note field
        eval(['hdm.',className,'.note=char([]);']);
    elseif comment
        eval(['hdm.',className,'.note=strvcat(hdm.',className,'.note,txtline);']);
    else
        [tag,txtline]=strtok(txtline,':');
        %  Special code for the case where a user has managed to get the first one or two
        %  characters of the coefficient set into a channel name.
        tag=deblank(strtok(tag,'-'));
        %  Remove leading ':' and leading and trailing white spaces
        txtline=deblank(txtline(length(txtline):-1:2));
        txtline=deblank(txtline(length(txtline):-1:1));
        
        %  Identify text fields
        if strcmp(className,'MultiSphere_Data') & strcmp(tag,'HEADPOS')
            %  This line has both numerical and text fields delimited by tabs,
            lasttab=max(strfind(txtline,char(9)));
            value=struct('HEADPOS',sscanf(txtline(1:lasttab),'%g'),...
                'NOMINAL',txtline((lasttab+1):length(txtline)));
            %headpos=sscanf(txtline(1:lasttab),'%g');
            %nominal=txtline((lasttab+1):length(txtline));
            %hdm.MultiSphere_Data.HEADPOS=struct('HEADPOS',headpos,'NOMINAL',nominal);
        elseif strcmp(className,'File_Info') | strcmp(className,'Model') | ...
                (strcmp(className,'MultiSphere_Data') & ...
                (strcmp(tag,'HEADSHAPE_FILE') | strcmp(tag,'SURFACE_TYPE')))
            % Remove leading and trailing blanks and the ':' terminator of the tag name.
            value=txtline;
        else
            value=sscanf(txtline,'%g');
        end
        eval(['hdm.',className,'.',tag,'=value;']);
    end
end
fclose(fid);
if classOpen
    fprintf(['readCTFhdm : End-of file with class open?\n',...
            '\t\tThere must be an error in %s\n'],fileName);
end
clear txtline tag value classOpen className comment fid;

if isfield(hdm,'MultiSphere_Data');
    %  Convert to a more convenient way to interpret info.
    SQUIDname=char(fieldnames(hdm.MultiSphere_Data));
    %  v6.0 .hdm files define SURFACE_TYPE and HEADPOS
    nonSQUIDname=strvcat('note','SEARCH_RADIUS','HEADSHAPE_FILE','SURFACE_TYPE','HEADPOS');
    
    %  Remove tags referring to something else
    for k=1:size(nonSQUIDname,1)
        index=strmatch(deblank(nonSQUIDname(k,:)),SQUIDname,'exact');
        if ~isempty(index)
            SQUIDname=SQUIDname(setdiff(1:length(SQUIDname),index),:);
        end
    end
    clear k index;
    
    SQUIDname=deblank(SQUIDname);
    sphereOrigin=zeros(3,size(SQUIDname,1));
    radius=zeros(1,size(SQUIDname,1));
    
    for k=size(SQUIDname,1):-1:1
        buff=getfield(hdm.MultiSphere_Data,deblank(SQUIDname(k,:)));
        sphereOrigin(:,k)=buff(1:3);
        radius(k)=buff(4);
        hdm.MultiSphere_Data=rmfield(hdm.MultiSphere_Data,deblank(SQUIDname(k,:)));
    end
    clear k buff;
    
    hdm.MultiSphere_Data.SQUIDname=SQUIDname;
    hdm.MultiSphere_Data.sphereOrigin=sphereOrigin;
    hdm.MultiSphere_Data.radius=radius;
end

% If there is no head-shape information in the file, print a warning.
if ~isfield(hdm,'Model') | (~isfield(hdm,'MultiSphere_Data') & ~isfield(hdm,'MEG_Sphere'))
    fprintf(['\nreadCTFhdm:  Structure hdm is missing fields Model, MEG_sphere, or ',...
            'MultiSphere_Data.\n\t\t\t Is file %s a head model file?\n\n'],fileName);
    missingField=1;
else
    missingField=0;
end

if isempty(printWarning)
    fprintf(['\nreadCTFhdm: You are reading a CTF head model for use with a software-application\n',...
            '\ttool that is not manufactured by VSM MedTech Ltd. and has not received marketing\n',...
            '\tclearance for clinical applications.  If CTF MEG data are processed by this tool,\n',...
            '\tthey should not be later employed for clinical and/or diagnostic purposes.\n']);
    if missingField==0;printWarning=1;end
end

return