function v = ft_version(cmd)

% FT_VERSION provides functionality for displaying version information on
% the current version of FieldTrip being used, as well us for updating the
% FieldTrip installation.
%
% To understand the different options this function provides, a little bit
% of background knowledge is needed about how FieldTrip deals with
% versions. FieldTrip does not use a fixed release system, but instead
% the 'release' version available from the FTP server is updated daily,
% through an automatic script. This release version is based on the
% development version used internally by the developers, which is updated
% more frequently. Development updates are tracked using Subversion (SVN),
% a system that assigns a revision number to each change to the code. The
% release version of FieldTrip contains a signature file (signature.md5),
% that contains, for each file in the release, the latest revision number
% for that file, along with an MD5 hash of the file's contents. A similar
% signature file is located on the FieldTrip webserver, which is updated
% whenever a developer makes a change to the codebase.
%
% This function uses both the local signature file (dating from whenever
% you downloaded FieldTrip) and the remote signature file (reflecting the
% latest version available--this updates more frequently than the daily FTP
% release) to determine (1) the version of FieldTrip you are using, (2)
% whether you have made local changes to your FieldTrip installation (by
% comparing the actual MD5 hashes of files to the ones stored in your local
% signature file), (3) whether there are new updates available (by
% comparing the local and remote signature files).
%
% Usage:
%   ft_version info             - display the latest revision number
%                                 present in your installation
%   ft_version full             - display full revision information about
%                                 your installation; this includes checking
%                                 for local and/or remote changes
%   ft_version update           - also displays full revision information;
%                                 additionally provides the option of
%                                 downloading new or updated files from the
%                                 code repository

% 'Hidden' options:
%   ft_version signature        - generate a new signature file, without
%                                 revision information about the files
%   ft_version fullsignature    - generate a new signature file, including
%                                 revision information; this can only be
%                                 used when you are using an SVN
%                                 installation of FieldTrip

% Copyright (C) 2012, Eelke Spaak
%
% This file is part of FieldTrip, see http://www.ru.nl/donders/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin<1
  cmd = 'info';
end

[ftpath, ~] = fileparts(mfilename('fullpath'));
signaturefile = fullfile(ftpath, 'signature.md5');
remotesignature = 'http://fieldtrip.fcdonders.nl/signature.md5';
repository = 'http://fieldtrip.googlecode.com/svn/trunk/';

% are we dealing with an SVN working copy of fieldtrip?
issvn = isdir(fullfile(ftpath, '.svn'));

switch cmd

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case 'info'
    % show the latest revision present in this copy of fieldtrip
    
    if issvn
      % use svn system call to determine latest revision
      olddir = pwd();
      cd(ftpath);
      [status,output] = system('svn info');
      cd(olddir);
      if status > 0
        error('you seem to have an SVN working copy of FieldTrip, yet ''svn info'' does not work as expected');
      end
      
      rev = regexp(output, 'Revision: (.*)', 'tokens', 'dotexceptnewline');
      rev = rev{1}{1};
      
    else
      % determine latest revision from the signature file
      fp = fopen(signaturefile, 'r');
      line = fgetl(fp); % just get first line, file should be ordered newest-first
      fclose(fp);
      
      rev = regexp(line, '[^\t]*\t([^\t])*\t.*', 'tokens');
      rev = rev{1}{1};
      
    end
    
    if nargout > 0
      v = rev;
    elseif issvn
      fprintf('\nThis is FieldTrip, version r%s (svn).\n\n', rev);
    else
      fprintf('\nThis is FieldTrip, version r%s.\n\n', rev);
    end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  case 'full'
    % show the latest revision present in this version of FieldTrip
    % and also check for possible remote and local changes
    
    if issvn
      % use svn status
      olddir = pwd();
      cd(ftpath);
      [status,output] = system('svn status');
      cd(olddir);
      if status > 0
        error('you seem to have an SVN working copy of FieldTrip, yet ''svn status'' does not work as expected');
      end
      
      if numel(output) > 1
        modFlag = 'M';
      else
        modFlag = 'U';
      end
      fprintf('\nThis is FieldTrip, version r%s%s (svn).\n\n', ft_version(), modFlag);
      
      % simply print the output of the svn status command
      if numel(output) > 1
        fprintf('This is the SVN status of your working copy:\n');
        fprintf(output);
        fprintf('\n');
      end
      
    else
      
      % fetch remote hash table
      str = urlread(remotesignature);
      [remHashes, remFiles, remRevs] = parseHashTable(str);

      % fetch local hash table
      str = fileread(signaturefile);
      [locHashes, locFiles, locRevs] = parseHashTable(str);

      % determine changes
      [remoteDeleted, remoteNew, remoteChanges,...
        localDeleted, localChanges] = findChanges(locHashes, locFiles,...
        remHashes, remFiles, ftpath);

      % print detailed version information
      if any(localChanges)
        modFlag = 'M'; % M for modified
      else
        modFlag = 'U'; % U for unmodified
      end
      fprintf('This is FieldTrip, version r%s%s.\n', ft_version(), modFlag);

      if any(localChanges)
        fprintf('\nthe following files have been locally modified:\n');
        inds = find(localChanges);
        for k = 1:numel(inds)
          fprintf('           r%sm    %s\n', locRevs{inds(k)}, locFiles{inds(k)});
        end
      end

      if any(localDeleted)
        fprintf('\nthe following files have been locally deleted:\n');
        inds = find(localDeleted);
        for k = 1:numel(inds)
          fprintf('           r%sd    %s\n', locRevs{inds(k)}, locFiles{inds(k)});
        end
      end

      if any(remoteChanges)
        fprintf('\na new version is available for the following files:\n');
        inds = find(remoteChanges);
        for k = 1:numel(inds)
          locInd = strcmp(locFiles, remFiles{inds(k)});
          fprintf('   r%s-->r%s     %s\n',...
            locRevs{locInd}, remRevs{inds(k)}, remFiles{inds(k)});
        end
      end

      if any(remoteDeleted)
        fprintf('\nthe following files are no longer part of FieldTrip:\n');
        files = locFiles(remoteDeleted);
        for k = 1:numel(files)
          fprintf('                     %s\n', files{k});
        end
      end

      if any(remoteNew)
        fprintf('\nthe following new files are available:\n');
        files = remFiles(remoteNew);
        revs = remRevs(remoteNew);
        for k = 1:numel(files)
          fprintf('           r%sa    %s\n', revs{k}, files{k});
        end
      end

      fprintf('\n');
      
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  case 'update'
    % get the remote changes, check that there are no conflicting local
    % changes and update the files
    
    if issvn
      
      fprintf('\nThis is an SVN working copy of FieldTrip, invoking ''svn update''...\n');
      olddir = pwd();
      cd(ftpath);
      system('svn update');
      cd(olddir);
      fprintf('\n');
      
    else
    
      % fetch remote hash table
      str = urlread(remotesignature);
      [remHashes, remFiles, remRevs] = parseHashTable(str);

      % fetch local hash table
      str = fileread(signaturefile);
      [locHashes, locFiles, locRevs] = parseHashTable(str);
      
      % determine changes
      [remoteDeleted, remoteNew, remoteChanges,...
        localDeleted, localChanges] = findChanges(locHashes, locFiles,...
        remHashes, remFiles, ftpath);
      
      % index all possible deletions and additions
      locFilesDelete = [];
      remFilesDownload = [];
      
      olddir = pwd();
      cd(ftpath);
      
      
      % keep track of whether local additions will possibly be overwritten
      overwriteFlag = 0;
      
      % first display information on what would be changed
      if any(remoteNew)
        fprintf('\nthe following new files will be added:\n');
        inds = find(remoteNew);
        for k = 1:numel(inds)
          if exist(remFiles{inds(k)}, 'file')
            msg = '!';
            overwriteFlag = 1;
          else
            msg = ' ';
          end
          
          % give this change a number and remember it
          remFilesDownload(end+1) = inds(k);
          locFilesDelete(end+1) = nan;
          
          fprintf('%s %3d.          r%s  %s\n', msg, numel(locFilesDelete), remRevs{inds(k)}, remFiles{inds(k)});
        end
      end
      
      if any(remoteDeleted)
        fprintf('\nthe following files will be deleted:\n');
        inds = find(remoteDeleted);
        for k = 1:numel(inds)
          if localChanges(inds(k))
            msg = '!';
            overwriteFlag = 1;
          else
            msg = ' ';
          end
          
          % give this change a number and remember it
          remFilesDownload(end+1) = nan;
          locFilesDelete(end+1) = inds(k);
          
          fprintf('%s %3d.                 %s\n', msg, numel(locFilesDelete), locFiles{inds(k)});
        end
      end
      
      if any(remoteChanges)
        fprintf('\nthe following files will be updated:\n');
        inds = find(remoteChanges);
        for k = 1:numel(inds)
          locInd = strcmp(locFiles, remFiles{inds(k)});
          if localChanges(locInd)
            msg = '!';
            overwriteFlag = 1;
          else
            msg = ' ';
          end
          
          % give this change a number and remember it
          remFilesDownload(end+1) = inds(k);
          locFilesDelete(end+1) = nan;
          
          fprintf('%s %3d.  r%s-->r%s  %s\n',...
            msg, numel(locFilesDelete), locRevs{locInd},...
            remRevs{inds(k)}, remFiles{inds(k)});
        end
      end
      
      assert(numel(remFilesDownload) == numel(locFilesDelete));
      
      if overwriteFlag
        fprintf(['\n!!NOTE: if you choose to include changes marked with a leading !,\n'...
          '        they will overwrite local changes or additions!\n']);
      end
      
      % prompt the user which changes should be incorporated
      while (true) % use while (true) as poor man's do{}while()
        changesToInclude = input(['\nwhich changes do you want to include?\n'...
          '([vector] for specific changes, ''all'', or [] for none)\n? '],'s');

        % check integrity of entered data
        if strcmp(changesToInclude, 'all')
          changesToInclude = 1:numel(remFilesDownload);
        else
          [changesToInclude,status] = str2num(changesToInclude);
        end
        
        if status && isnumeric(changesToInclude) && all(changesToInclude > 0)...
            && (all(changesToInclude <= numel(remFilesDownload)) || ...
              changesToInclude == Inf)
            changesToInclude = unique(changesToInclude);
            break;
        end
      end
      
      % include only those changes requested in the actual update
      locFilesDelete = locFilesDelete(changesToInclude);
      locFilesDelete(isnan(locFilesDelete)) = [];
      remFilesDownload = remFilesDownload(changesToInclude);
      remFilesDownload(isnan(remFilesDownload)) = [];
      
      % delete all local files that should be deleted
      for k = 1:numel(locFilesDelete)
        ind = locFilesDelete(k);
        fprintf('deleting file %s...\n', locFiles{ind});
        delete(locFiles{ind});
        % also remove entry from the hash table cell arrays
        locFiles(ind) = [];
        locHashes(ind) = [];
        locRevs(ind) = [];
      end
      
      % add new files and/or update changed files
      for k = 1:numel(remFilesDownload)
        ind = remFilesDownload(k);
        fprintf('downloading file %s...\n', remFiles{ind});
        newFile = urlread([repository remFiles{ind}]);
        fp = fopen(remFiles{ind}, 'w');
        fwrite(fp, newFile);
        fclose(fp);
        
        % update local hash table
        locInd = strcmp(locFiles, remFiles{ind});
        if any(locInd)
          % file was already present, update hash
          locHashes{locInd} = remHashes{ind};
          locRevs{locInd} = remRevs{ind};
        else
          % new file, add hash
          locFiles{end} = remFiles{ind};
          locHashes{end} = remHashes{ind};
          locRevs{end} = remRevs{ind};
        end
        
      end
      
      cd(olddir);
      
      if ~isempty(changesToInclude)
        % output new hash table
        fp = fopen(signaturefile,'w');
        outputHashTable(fp, locHashes, locFiles, locRevs);
        fclose(fp);

        fprintf('update finished.\n');
      else
        fprintf('nothing updated.\n');
      end
       
      ft_version();
      
    end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  case 'signature'
    
    if exist(signaturefile, 'file')
      error(sprintf([signaturefile ' already exists.\nIf you are sure you want to re-create '...
        'the signature file, you will have to manually delete it.\nNote that this means '...
        'that local edits cannot be tracked anymore, and might be overwritten without notice!']));
    end
    
    fprintf('\nwriting signature file to signature-local.md5 (this might take a few minutes)...');
    hashFile = fopen(signaturefile, 'w');
    hashAll(hashFile, ftpath, ftpath);
    fclose(hashFile);
    fprintf('done.\n\n');
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  case 'fullsignature'
    
    if exist(signaturefile, 'file')
      error(sprintf([signaturefile ' already exists.\nIf you are sure you want to re-create '...
        'the signature file, you will have to manually delete it.\nNote that this means '...
        'that local edits cannot be tracked anymore, and might be overwritten without notice!']));
    end
    
    if ~issvn
      error('you can only generate a full signature with revision info if you are using an SVN working copy of FieldTrip');
    end

    fprintf('\nwriting signature file to signature-local.md5 (this might take a few minutes)...');
    hashFile = fopen(fullfile(ftpath, 'signature-local.md5'), 'w');
    hashAll(hashFile, ftpath, ftpath, 1);
    fclose(hashFile);
    fprintf('done.\n\n');

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% locHashes and locFiles represent the hashes and filenames of entries in
% the local hash table, respectively. remHashes and remFiles represent the
% hashes and filenames of entries in the remote hash table, respectively.
% ftpath is the FT root path. Returned are logical index vectors into the
% remFiles cell array (for remoteChanges and remoteNew) or into the
% locFiles cell array (rest).
function [remoteDeleted, remoteNew, remoteChanges,...
  localDeleted, localChanges] = findChanges(locHashes, locFiles,...
  remHashes, remFiles, ftpath)

  fprintf('\ncomparing local and remote (latest) hashes, this might take a minute...\n');

  remoteDeleted = false(size(locFiles));
  localChanges = false(size(locFiles));
  remoteChanges = false(size(remFiles));
  localDeleted = false(size(locFiles));
  % note: no keeping track of local additions (not desirable probably)

  % compare them by looping over all files in local table
  for k = 1:numel(locFiles)
    ind = strcmp(remFiles, locFiles{k});

    if sum(ind) > 1
      % this should never happen; means remote table is corrupt
      warning('more than one entry found in remote hash table for file %s', locFiles{k});
    end

    if ~any(ind)
      % file not found in remote table, should be deleted locally
      remoteDeleted(k) = 1;
      continue;
    end

    if ~strcmp(remHashes{ind}, locHashes{k})
      % hash remote different from hash in local table
      remoteChanges(ind) = 1;
    end

    if ~exist(fullfile(ftpath, locFiles{k}), 'file')
      localDeleted(k) = 1;
      continue;
    end

    if (nargout > 4) % only check for local changes if requested, slow step
      if ~strcmp(CalcMD5(fullfile(ftpath, locFiles{k}), 'File'), locHashes{k})
        % hash of actual file different from hash in local table
        localChanges(k) = 1;
      end
    end
  end

  % only remote additions have not yet been determined, do so now
  remoteNew = false(size(remFiles));
  [~,inds] = setdiff(remFiles, locFiles);
  remoteNew(inds) = 1;
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If str is a tab-separated string, with records separated by newlines,
% parseHashTable(str) will return the three columns in the file: the first
% column will be returned as hashes; the second as revs, the third as
% files.
function [hashes,files,revs] = parseHashTable(str)

  lines = regexp(str, '\n', 'split');
  if isempty(lines{end})
    lines = lines(1:end-1);
  end

  hashes = cell(1, numel(lines));
  files = cell(1, numel(lines));
  revs = cell(1, numel(lines));

  for k = 1:numel(lines)
    tmp = regexp(lines{k}, '\t', 'split');
    hashes{k} = tmp{1};
    revs{k} = tmp{2};
    files{k} = tmp{3};
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hashAll(filePointer, path, basepath, withRev)

  if nargin < 4
    withRev = 0;
  end
  
  % these files are always excluded from the update and hashing routines
  file_excludes = {'.', '..', '.svn', 'test', 'signature.md5', 'signature-local.md5'};

  list = dir(path);
  for k = 1:numel(list)
    if ~any(strcmp(list(k).name, file_excludes))
      filename = [path '/' list(k).name];
      if (list(k).isdir)
        hashAll(filePointer, filename, basepath, withRev);
      else
        md5 = CalcMD5(filename,'File');
        
        if withRev
          [status, output] = system(['svn info ' filename]);
          rev = regexp(output, 'Last Changed Rev: (.*)', 'tokens', 'dotexceptnewline');
          if ~isempty(rev) && ~isempty(rev{1})
            rev = rev{1}{1};
          else
            rev = 'UNKNOWN';
          end
        else
          rev = 'UNKNOWN';
        end
        filename = filename((numel(basepath)+2):end); % strip off the base path
        fprintf(filePointer, '%s\t%s\t%s\n',md5,rev,filename);
      end
    end
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outputHashTable(filePointer, hashes, files, revs)

  % sort the entries (highest revs should always be at the top)
  numRevs = cellfun(@str2double, revs);
  [sortedRevs,inds] = sort(numRevs, 'descend');
  lastnan = find(isnan(sortedRevs), 1, 'last'); % put nans (=UNKNOWN) at the end
  inds = [inds(lastnan+1:end) inds(1:lastnan)];

  hashes = hashes(inds);
  files = files(inds);
  revs = revs(inds);

  for k = 1:numel(hashes)
    fprintf(filePointer, '%s\t%s\t%s\n', hashes{k}, revs{k}, files{k});
  end

end