function spm_check_installation(action)
% Check SPM installation
% FORMAT spm_check_installation('basic')
% Perform a superficial check of SPM installation [default].
%
% FORMAT spm_check_installation('full')
% Perform an in-depth diagnostic of SPM installation.
%
% FORMAT spm_check_installation('build')
% Build signature of SPM distribution as used by previous option.
% (for developers)
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id$

if isdeployed, return; end

%-Select action
%--------------------------------------------------------------------------
if ~nargin, action = 'basic'; end
switch action
    case 'basic'
        check_basic;
    case 'full'
        check_full;
    case 'build'
        build_signature;
    otherwise
        error('Unknown action to perform.');
end

%==========================================================================
% FUNCTION check_basic
%==========================================================================
function check_basic

%-Minimal MATLAB version required
%--------------------------------------------------------------------------
try
    v = spm_matlab_version_chk('7.1');
catch
    error('Where is spm_matlab_version_chk.m?');
end
if v < 0
    error([...
        'SPM8 requires MATLAB 7.1 onwards in order to run.\n'...
        'This MATLAB version is %s\n'], version);
end

%-Check installation
%--------------------------------------------------------------------------
spm('Ver','',1);
d = spm('Dir');

%-Check the search path
%--------------------------------------------------------------------------
p = textscan(path,'%s','delimiter',pathsep); p = p{1};
if ~ismember(lower(d),lower(p))
    error(sprintf([...
        'You do not appear to have the MATLAB search path set up\n'...
        'to include your SPM8 distribution. This means that you\n'...
        'can start SPM in this directory, but if your change to\n'...
        'another directory then MATLAB will be unable to find the\n'...
        'SPM functions. You can use the editpath command in MATLAB\n'...
        'to set it up.\n'...
        '    addpath %s\n'...
        'For more information, try typing the following:\n'...
        '    help path\n    help editpath'],d));
end

%-Ensure that the original release - as well as the updates - was installed
%--------------------------------------------------------------------------
if ~exist(fullfile(d,'spm_Npdf.m'),'file') % File that should not have changed
    if isunix
        error(sprintf([...
            'There appears to be some problem with the installation.\n'...
            'The original spm8.zip distribution should be installed\n'...
            'and the updates installed on top of this. Unix commands\n'...
            'to do this are:\n'...
            '   unzip spm8.zip\n'...
            '   unzip -o spm8_updates_r????.zip -d spm8\n'...
            'You may need help from your local network administrator.']));
    else
        error(sprintf([...
            'There appears to be some problem with the installation.\n'...
            'The original spm8.zip distribution should be installed\n'...
            'and the updates installed on top of this. If in doubt,\n'...
            'consult your local network administrator.']));
    end
end

%-Ensure that the files were unpacked correctly
%--------------------------------------------------------------------------
if ispc
    try
        t = load(fullfile(d,'Split.mat'));
    catch
        error(sprintf([...
            'There appears to be some problem reading the MATLAB .mat\n'...
            'files from the SPM distribution. This is probably\n'...
            'something to do with the way that the distribution was\n'...
            'unpacked. If you used WinZip, then ensure that\n'...
            'TAR file smart CR/LF conversion is disabled\n'...
            '(under the Miscellaneous Configuration Options).']));
    end
    if ~exist(fullfile(d,'toolbox','DARTEL','diffeo3d.c'),'file'),
        error(sprintf([...
            'There appears to be some problem with the installation.\n'...
            'This is probably something to do with the way that the\n'...
            'distribution was unbundled from the original .zip files.\n'...
            'Please ensure that the files are unpacked so that the\n'...
            'directory structure is retained.']));
    end
end

%-Check the MEX files
%--------------------------------------------------------------------------
try
    feval(@spm_bsplinc,1,ones(1,6));
catch
    error([...
        'SPM uses a number of MEX files, which are compiled functions.\n'...
        'These need to be compiled for the various platforms on which SPM\n'...
        'is run. At the FIL, where SPM is developed, the number of\n'...
        'computer platforms is limited.  It is therefore not possible to\n'...
        'release a version of SPM that will run on all computers. See\n'...
        '   %s%csrc%cMakefile and\n'...
        '   http://en.wikibooks.org/wiki/SPM#Installation\n'...
        'for information about how to compile mex files for %s\n'...
        'in MATLAB %s.'],...
        d,filesep,filesep,computer,version);
end

%==========================================================================
% FUNCTION check_full
%==========================================================================
function check_full

%-Say Hello
%--------------------------------------------------------------------------
disp( ' ___  ____  __  __                                           ' );
disp( '/ __)(  _ \(  \/  )                                          ' );
disp( '\__ \ )___/ )    (   Statistical Parametric Mapping          ' );
disp(['(___/(__)  (_/\/\_)  SPM - http://www.fil.ion.ucl.ac.uk/spm/ ']);
fprintf('\n');

%-Detect SPM directory
%--------------------------------------------------------------------------
SPMdir = which('spm.m','-ALL');
if isempty(SPMdir)
    fprintf('SPM is not in your MATLAB path.\n');
    return;
elseif numel(SPMdir) > 1
    fprintf('SPM seems to appear in several different folders:\n');
    for i=1:numel(SPMdir)
        fprintf('  * %s\n',SPMdir{i});
    end
    fprintf('Remove all but one with ''editpath'' or ''spm_rmpath''.\n');
    return;
else
    fprintf('SPM is installed in: %s\n',fileparts(SPMdir{1}));
end
SPMdir = fileparts(SPMdir{1});

%-Detect SPM version and revision number
%--------------------------------------------------------------------------
v = struct('Name','','Version','','Release','','Date','');
try
    fid = fopen(fullfile(SPMdir,'Contents.m'),'rt');
    if fid == -1
        fprintf('Cannot open ''%s'' for reading.\n',fullfile(SPMdir,'Contents.m'));
        return;
    end
    l1 = fgetl(fid); l2 = fgetl(fid);
    fclose(fid);
    l1 = strtrim(l1(2:end)); l2 = strtrim(l2(2:end));
    t  = textscan(l2,'%s','delimiter',' '); t = t{1};
    v.Name = l1; v.Date = t{4};
    v.Version = t{2}; v.Release = t{3}(2:end-1);
catch
    fprintf('Cannot obtain SPM version & revision number.\n');
    return;
end
fprintf('SPM version is %s (%s, %s)\n', ...
    v.Release,v.Version,strrep(v.Date,'-',' '));

%-Detect SPM toolboxes
%--------------------------------------------------------------------------
officials = {'Beamforming', 'DARTEL', 'dcm_meeg', 'DEM', 'FieldMap', ...
    'HDW', 'MEEGtools', 'mixture', 'Neural_Models', 'Seg', 'spectral', ...
    'SRender'};
dd = dir(fullfile(SPMdir,'toolbox'));
dd = {dd([dd.isdir]).name};
dd(strmatch('.',dd)) = [];
dd = setdiff(dd,officials);
fprintf('SPM toolboxes:');
for i=1:length(dd)
    fprintf(' %s',dd{i});
end
if isempty(dd), fprintf(' none'); end
fprintf('\n');

%-Detect MATLAB & toolboxes
%--------------------------------------------------------------------------
fprintf('MATLAB is installed in: %s\n',matlabroot);
fprintf('MATLAB version is %s\n',version);
fprintf('MATLAB toolboxes: '); hastbx = false;
if license('test','signal_toolbox') && ~isempty(ver('signal'))
    vtbx = ver('signal'); hastbx = true;
    fprintf('signal (v%s) ',vtbx.Version);
end
if license('test','image_toolbox') && ~isempty(ver('images'))
    vtbx = ver('images'); hastbx = true;
    fprintf('images (v%s) ',vtbx.Version);
end
if license('test','statistics_toolbox') && ~isempty(ver('stats'))
    vtbx = ver('stats'); hastbx = true;
    fprintf('stats (v%s)',vtbx.Version);
end
if ~hastbx, fprintf('none.'); end
fprintf('\n');

%-Detect Platform and Operating System
%--------------------------------------------------------------------------
[C, maxsize] = computer;
fprintf('Platform: %s (maxsize=%d)\n', C, maxsize);
if ispc
   platform = [system_dependent('getos'),' ',system_dependent('getwinsys')];
elseif exist('ismac') && ismac
    [fail, input] = unix('sw_vers');
    if ~fail
    platform = strrep(input, 'ProductName:', '');
    platform = strrep(platform, sprintf('\t'), '');
    platform = strrep(platform, sprintf('\n'), ' ');
    platform = strrep(platform, 'ProductVersion:', ' Version: ');
    platform = strrep(platform, 'BuildVersion:', 'Build: ');
    else
        platform = system_dependent('getos');
    end
else    
   platform = system_dependent('getos');
end
fprintf('OS: %s\n', platform);

%-Detect Java
%--------------------------------------------------------------------------
fprintf('%s\n', version('-java'));
fprintf('Java support: ');
level = {'jvm', 'awt', 'swing', 'desktop'};
for i=1:numel(level)
    if isempty(javachk(level{i})), fprintf('%s ',level{i}); end
end
fprintf('\n');

%-Detect Monitor(s)
%--------------------------------------------------------------------------
M = get(0,'MonitorPositions');
fprintf('Monitor(s):');
for i=1:size(M,1)
    fprintf(' [%d %d %d %d]',M(i,:));
end
fprintf(' (%dbit)\n', get(0,'ScreenDepth'));

%-Detect OpenGL rendering
%--------------------------------------------------------------------------
S =  opengl('data');
fprintf('OpenGL version: %s',S.Version);
if S.Software, fprintf('(Software)\n'); else fprintf('(Hardware)\n'); end
fprintf('OpenGL renderer: %s (%s)\n',S.Vendor,S.Renderer);

%-Detect MEX setup
%--------------------------------------------------------------------------
fprintf('MEX extension: %s\n',mexext);
try
    cc = mex.getCompilerConfigurations('C','Selected');
    if ~isempty(cc)
        cc = cc(1); % can be C or C++
        fprintf('C Compiler: %s (%s).\n', cc.Name, cc.Version);
        fprintf('C Compiler settings: %s (''%s'')\n', ...
            cc.Details.CompilerExecutable, cc.Details.OptimizationFlags);
    else
        fprintf('No C compiler is selected (see mex -setup)\n');
    end
end
try
    [sts, m] = fileattrib(fullfile(SPMdir,'src'));
    m = [m.UserRead m.UserWrite m.UserExecute ...
         m.GroupRead m.GroupWrite m.GroupExecute ...
         m.OtherRead m.OtherWrite m.OtherExecute];
    r = 'rwxrwxrwx'; r(~m) = '-';
    fprintf('C Source code permissions: dir %s, ', r);
    [sts, m] = fileattrib(fullfile(SPMdir,'src','spm_resels_vol.c'));
    m = [m.UserRead m.UserWrite m.UserExecute ...
         m.GroupRead m.GroupWrite m.GroupExecute ...
         m.OtherRead m.OtherWrite m.OtherExecute];
    r = 'rwxrwxrwx'; r(~m) = '-';
    fprintf('file %s\n',r);
end

%-Get file details for local SPM installation
%--------------------------------------------------------------------------
fprintf('%s\n',repmat('-',1,70));
l = generate_listing(SPMdir);
fprintf('%s %40s\n','Parsing local installation...','...done');

%-Get file details for most recent public version
%--------------------------------------------------------------------------
fprintf('Downloading SPM information...');
url = sprintf('http://www.fil.ion.ucl.ac.uk/spm/software/%s/%s.xml',...
    lower(v.Release),lower(v.Release));
try
    p = [tempname '.xml'];
    urlwrite(url,p);
catch
    fprintf('\nCannot access URL %s\n',url);
    return;
end
fprintf('%40s\n','...done');

%-Parse it into a Matlab structure
%--------------------------------------------------------------------------
fprintf('Parsing SPM information...');
tree = xmlread(p);
delete(p);
r    = struct([]); ind = 1;
if tree.hasChildNodes
    for i=0:tree.getChildNodes.getLength-1
        m = tree.getChildNodes.item(i);
        if strcmp(m.getNodeName,'signature')
            m = m.getChildNodes;
            for j=0:m.getLength-1
                if strcmp(m.item(j).getNodeName,'file')
                    n = m.item(j).getChildNodes;
                    for k=0:n.getLength-1
                        o = n.item(k);
                        if ~strcmp(o.getNodeName,'#text')
                            try
                                s = char(o.getChildNodes.item(0).getData);
                            catch
                                s = '';
                            end
                            switch char(o.getNodeName)
                                case 'name'
                                    r(ind).file = s;
                                case 'id'
                                    r(ind).id = str2num(s);
                                otherwise
                                    r(ind).(char(o.getNodeName)) = s;
                            end
                        end
                    end
                    ind = ind + 1;
                end
            end
        end
    end
end
fprintf('%44s\n','...done');

%-Compare local and public versions
%--------------------------------------------------------------------------
fprintf('%s\n',repmat('-',1,70));
compare_versions(l,r);
fprintf('%s\n',repmat('-',1,70));

%==========================================================================
% FUNCTION compare_versions
%==========================================================================
function compare_versions(l,r)

%-Uniformise file names
%--------------------------------------------------------------------------
a = {r.file}; a = strrep(a,'\','/'); [r.file] = a{:};
a = {l.file}; a = strrep(a,'\','/'); [l.file] = a{:};

%-Look for missing or unknown files
%--------------------------------------------------------------------------
[x,ir,il] = setxor({r.file},{l.file});
if isempty([ir il])
    fprintf('No missing or unknown files\n');
else
    if ~isempty(ir)
        fprintf('File(s) missing in your installation:\n');
    end
    for i=1:length(ir)
        fprintf(' * %s\n', r(ir(i)).file);
    end
    if ~isempty(ir) && ~isempty(il), fprintf('%s\n',repmat('-',1,70)); end
    if ~isempty(il)
        fprintf('File(s) not part of the current version of SPM:\n');
    end
    for i=1:length(il)
        fprintf(' * %s\n', l(il(i)).file);
    end
end

fprintf('%s\n',repmat('-',1,70));

%-Look for local changes or out-of-date files
%--------------------------------------------------------------------------
[tf, ir] = ismember({r.file},{l.file});
dispc = true;
for i=1:numel(r)
    if tf(i)
        if ~isempty(r(i).id) && (isempty(l(ir(i)).id) || (r(i).id ~= l(ir(i)).id))
            dispc = false;
            fprintf('File %s is not up to date (r%d vs r%d)\n', ...
                r(i).file, r(i).id, l(ir(i)).id);
        end
        if ~isempty(r(i).md5) && ~isempty(l(ir(i)).md5) && ~strcmp(r(i).md5,l(ir(i)).md5)
            dispc = false;
            fprintf('File %s has been edited (checksum mismatch)\n', r(i).file);
        end
    end
end
if dispc
    fprintf('No local change or out-of-date files\n');
end

%==========================================================================
% FUNCTION build_signature
%==========================================================================
function build_signature

[v,r] = spm('Ver','',1);
d     = spm('Dir');
l     = generate_listing(d);

fprintf('Saving %s signature in %s\n',...
    v, fullfile(pwd,sprintf('%s_signature.xml',v)));
fid = fopen(fullfile(pwd,sprintf('%s_signature.xml',v)),'wt');
fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<!-- Signature for %s r%s -->\n',v,r);
fprintf(fid,'<signature>\n');
for i=1:numel(l)
    fprintf(fid,'  <file>\n');
    fprintf(fid,'    <name>%s</name>\n',l(i).file);
    fprintf(fid,'    <id>%d</id>\n',l(i).id);
    fprintf(fid,'    <date>%s</date>\n',l(i).date);
    fprintf(fid,'    <md5>%s</md5>\n',l(i).md5);
    fprintf(fid,'  </file>\n');
end
fprintf(fid,'</signature>\n');
fclose(fid);

%==========================================================================
% FUNCTION generate_listing
%==========================================================================
function l = generate_listing(d,r,l)

if nargin < 2
    r = '';
end
if nargin < 3
    l = struct([]);
end

%-List content of folder
%--------------------------------------------------------------------------
ccd = fullfile(d,r);
fprintf('%-10s: %58s','Directory',ccd(max(1,length(ccd)-57):end));
dd  = dir(ccd);
f   = {dd(~[dd.isdir]).name};
dispw = false;

for i=1:length(f)
    [p,name,ext] = fileparts(f{i});
    info = struct('file','', 'id',[], 'date','', 'md5','');
    if ismember(ext,{'.m','.man','.txt','.xml','.c','.h',''})
        info = extract_info(fullfile(ccd,f{i}));
    end
    info.file = fullfile(r,f{i});
    try
        try
            info.md5 = md5sum(fullfile(ccd,f{i}));
        catch
            %[s, info.md5] = system(['md5sum "' fullfile(ccd,f{i}) '"']);
            %info.md5 = strtok(info.md5);
        end
    end
    if isempty(l), l = info;
    else l(end+1) = info;
    end
    if isempty(r) && strcmp(ext,'.m')
        w = which(f{i},'-ALL');
        if numel(w) > 1 && ~strcmpi(f{i},'Contents.m')
            if ~dispw, fprintf('\n'); end
            dispw = true;
            fprintf('File %s appears %d times in your MATLAB path:\n',f{i},numel(w));
            for j=1:numel(w)
                if j==1 && isempty(strmatch(d,w{1}))
                    fprintf('  %s (SHADOWING)\n',w{1});
                else
                    fprintf('  %s\n',w{j});
                end
            end
        end
    end
end
if ~dispw, fprintf('%s',repmat(sprintf('\b'),1,70)); end

%-Recursively extract subdirectories
%--------------------------------------------------------------------------
dd = {dd([dd.isdir]).name};
dd(strmatch('.',dd)) = [];
for i=1:length(dd)
    l = generate_listing(d,fullfile(r,dd{i}),l);
end

%==========================================================================
% FUNCTION extract_info
%==========================================================================
function svnprops = extract_info(f)
%Extract Subversion properties ($Id$ tag)

svnprops = struct('file',f, 'id',[], 'date','', 'md5','');

fp  = fopen(f,'rt');
str = fread(fp,Inf,'*uchar');
fclose(fp);
str = char(str(:)');

r = regexp(str,['\$Id: (?<file>\S+) (?<id>[0-9]+) (?<date>\S+) ' ...
                '(\S+Z) (?<author>\S+) \$'],'names');
                
if isempty(r)
    %fprintf('\n%s has no SVN Id.\n',f);
else
    svnprops.file = r(1).file;
    svnprops.id   = str2num(r(1).id);
    svnprops.date = r(1).date;
    [p,name,ext]  = fileparts(f);
    if ~strmatch(svnprops.file,[name ext],'exact')
        warning('SVN Id for %s does not match filename.',f);
    end
end
if numel(r) > 1
    %fprintf('\n%s has several SVN Ids.\n',f);
end
