function [condNumbers,condLabels] = read_ctf_cls(fname)

% READ_CTF_CLS reads the classification file from a CTF dataset

% Copyright (C) 2003, Ole Jensen
%
% $Log: read_ctf_cls.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.2  2006/04/20 12:07:31  roboos
% fixed bug when a class had no trials associated to it, changed output format for condNumbers from cell-array with cells into a cell-array with numeric arrays
%
% Revision 1.1  2004/09/27 15:34:30  roboos
% copy from the readClassFile that Ole wrote
% added help comments and cleaned up white space
% no code changes
%

condNumbers = [];
% condLabels = [];

fid = fopen(fname,'r');

if fid == -1
  condNumbers = [];
  condLabels = [];
  return
end

nCondition = 0;
readBad = 0;
readList = 0;
S2 = '*';
S1 = '*';
while ~isempty(S1)
  S3 = S2;
  S2 = S1;
  S1 =fscanf(fid,'%s',1);

  if readList
    if ~isempty(S1) & ~isempty(str2num(S1(2:end)))
      k = k + 1;
      condTmp = [condTmp 1+str2num(S1(2:end))];
    else
      readList = 0;
    end
    condNumbers{nCondition} = condTmp;
  end

  if strcmp(S2,'NAME:')
    % New condition found!
    % fprintf('%s\n',S1);
    nCondition = nCondition+1;
    condLabels(nCondition) = {S1} ;
  end

  if strcmp(S1,'NUMBER') & strcmp(S2,'TRIAL')
    if ~isempty(S1)
      readList = 1;
      k = 0;
      condTmp = [];
    else
      readList = 0;
    end
  end

end % ~isempty(S1)

