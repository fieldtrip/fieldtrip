function [condNumbers,condLabels] = read_ctf_cls(fname)

% READ_CTF_CLS reads the classification file from a CTF dataset

% Copyright (C) 2003, Ole Jensen
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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

