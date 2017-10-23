function test_bug2093

% WALLTIME 00:20:00
% MEM 3gb

% TEST ft_read_header ft_read_data ft_read_event read_nex_header read_nex_data read_plexon_nex

% Compares output of same test NEX file from ft_read_header,
% ft_read_data, and ft_read_event using old (master branch) and new
% (nexhandling branch) versions of the FieldTrip code

testnexfile = dccnpath('/home/common/matlab/fieldtrip/data/test/original/lfp/plexon/p213parall.nex');
testmatfile = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2093.mat');

%% first get output of old code
% if ~exist(testmatfile, 'file')
%   % run only once, this is to be done using FT version older than 9/9/2016,
%   % last commit before changes:  c8de94a3e13df7d10737381f1070a0d467c9bfbb
%   old.hdr = ft_read_header(testnexfile);
%   old.dat = ft_read_data(testnexfile);
%   old.evt = ft_read_event(testnexfile);
%   
%   save(testmatfile, 'old')
% else
  % all future test executions with newer FT code
  load(testmatfile)
% end

%% then get output of new code
% new.hdr = ft_read_header(testnexfile);
% new.dat = ft_read_data(testnexfile);
new.evt = ft_read_event(testnexfile);

%% then compare the 2
% make sure the old header and data are identical, since I haven't changed them yet
% assert(isequal(old.hdr,new.hdr))
% assert(isequaln(old.dat,new.dat)) % data padded with NaNs should still be equal

%% make sure all old events are still present, despite adding more
assert(isequal(fieldnames(old.evt), fieldnames(new.evt)))
fnames = fieldnames(old.evt);
oldcell = cell(1,numel(fnames));
newcell = cell(1,numel(fnames));
mtnew = false(numel(new.evt),numel(fnames));
for fieldn = 1:numel(fnames)
  if ischar(old.evt(1).(fnames{fieldn})) && ...
      ischar(new.evt(1).(fnames{fieldn}))
    oldcell{fieldn} = {old.evt.(fnames{fieldn})}';
    newcell{fieldn} = {new.evt.(fnames{fieldn})}';
  else
    oldcell{fieldn} = vertcat(old.evt.(fnames{fieldn}));
    newcell{fieldn} = vertcat(new.evt.(fnames{fieldn}));
  end
  assert(isequal(class(oldcell{fieldn}), class(newcell{fieldn})))
  
  if ~isequal(size(oldcell{fieldn},1), numel(old.evt))
    error('This strategy won''t work.')
  end
  if ~isequal(size(newcell{fieldn},1), numel(new.evt))
    for evn = 1:numel(new.evt)
      mtnew(evn,fieldn) = isempty(new.evt(evn).(fnames{fieldn}));
    end
    if ~isequal(size(newcell{fieldn},1), sum(~mtnew(:,fieldn)))
      error('Empty fields don''t account for difference in event #s.')
    end
  end
end

oldevttable = struct2table(old.evt);
% Empty values in numeric variables cause the column to be represented as a
% cell array when the structure is converted to a table, which is not
% supported by ismember.  The empty values should only be in the new event
% types I added support for, so they can be safely excluded from the table
% when verifying that all old event rows are contained in the new event
% data.
newrows = ~any(mtnew,2);
newevttable = table; %#ok<*AGROW>
for fieldn = 1:numel(fnames)
  if isequal(size(newcell{fieldn},1), numel(newrows))
    newevttable = [newevttable table(newcell{fieldn}(newrows,:), ...
      'VariableNames',fnames(fieldn))];
  elseif isequal(size(newcell{fieldn},1), sum(newrows))
    newevttable = [newevttable table(newcell{fieldn}, ...
      'VariableNames',fnames(fieldn))];
  else
    error('There may be different empty cells in different event fields.')
  end
end

assert(all(ismember(oldevttable,newevttable)))
disp('All old events found in new event output.')

