function index = documentationindex(filename)

fieldtripdefs

funname = {
  'definetrial'
  'artifact_eog'
  'artifact_jump'
  'artifact_muscle'
  'rejectartifact'
  'preprocessing'
  'rejectvisual'
  'appenddata'
  'resampledata'
  'channelrepair'
  'recodeevent'
  'redefinetrial'
  'read_header'
  'read_data'
  'read_event'
  'read_mri'
  'read_spike'
  'write_data'
  'write_event'

  'timelockanalysis'
  'timelockgrandaverage'

  'freqanalysis'
  'freqanalysis_mtmfft'
  'freqanalysis_mtmwelch'
  'freqanalysis_mtmconvol'
  'freqanalysis_wltconvol'
  'freqanalysis_tfr'
  'freqbaseline'
  'freqgrandaverage'
  'freqdescriptives'

  'dipolefitting'
  'dipolesimulation'
  'sourceanalysis'
  'sourcegrandaverage'
  'sourcedescriptives'
  'sourceinterpolate'
  'prepare_localspheres'
  'prepare_singleshell'
  'prepare_bemmodel'
  'prepare_leadfield'
  'prepare_atlas'
  'volumelookup'
  'volumenormalise'
  'volumesegment'

  'timelockstatistics'
  'freqstatistics'
  'sourcestatistics'

  'spikedownsample'
  'spikesplitting'

  'rt_asaproxy'
  'rt_brainampproxy'
  'rt_ctfproxy'
  'rt_fileproxy'
  'rt_neuralynxproxy'
  'rt_plexonproxy'
  'rt_signalproxy'

  'prepare_layout'
  'layoutplot'
  'topoplot'
  'singleplotER'
  'singleplotTFR'
  'topoplotER'
  'topoplotTFR'
  'multiplotER'
  'multiplotTFR'
  'sourceplot'
  };

ncfg  = 0;
index = {};

for j=1:length(funname)
  str = help(funname{j})
  str = tokenize(str, 10);

  % compact the help description
  for i=2:(length(str)-1)
    prevline = str{i-1};
    thisline = str{i  };
    nextline = str{i+1};
    if length(thisline)<5 || length(prevline)<5
      continue
    end
    try
      if ~isempty(regexp(prevline, '^ *cfg')) && isempty(regexp(thisline, '^ *cfg')) && ~all(thisline(1:3)==' ')
        % do not concatenate, this line starts a new paragraph
      elseif ~isempty(regexp(prevline, '^ *cfg')) && all(thisline(1:5)==' ')
        % concatenate a multiline cfg description
        thisline = cat(2, prevline, thisline);
        prevline = '';
      elseif ~all(prevline(1:5)==' ') && ~all(thisline(1:5)==' ') && isempty(regexp(thisline, '^ *cfg'))
        % concatenate the lines of a paragraph
        thisline = cat(2, prevline, thisline);
        prevline = '';
      elseif isempty(regexp(prevline, '^ *cfg')) && ~isempty(regexp(thisline, '^  cfg'))
       % previous line is a paragraph, this line starts with "cfg" but has no extra space in front of it
       % so assume that the cfg is part of the running text in the paragraph and conactenate the lines
       thisline = cat(2, prevline, thisline);
       prevline = '';
      end
    catch
      disp(lasterr);
      disp(thisline);
    end
    str{i-1} = prevline;
    str{i  } = thisline;
    str{i+1} = nextline;
  end
  for i=1:length(str)
    if length(str{i})>1
      % remove double spaces
      dum = findstr(str{i}, '  ');
      str{i}(dum) = [];
    end
    while ~isempty(str{i}) && str{i}(1)==' '
      % remove spaces at the begin of the line
      str{i}(1) = [];
    end
    while ~isempty(str{i}) && str{i}(end)==' '
      % remove spaces at the end of the line
      str{i}(1) = [];
    end
  end
  for i=1:length(str)
    if ~isempty(regexp(str{i}, '^ *cfg.[a-zA-Z0-9_\.]*'))
      ncfg = ncfg+1;
      index{ncfg,1} = funname{j};
      dum = regexp(str{i}, 'cfg.[a-zA-Z0-9_\.]*', 'match');
      index{ncfg,2} = dum{1};
      dum = str{i};
      while length(dum)>0 && dum(1)~=' '
        dum = dum(2:end);
      end
      while length(dum)>0 && (dum(1)=='=' || dum(1)==' ')
        dum = dum(2:end);
      end
      index{ncfg,3} = dum;
      dum1 = index{ncfg,1};
      dum1(end+1:30) = ' ';
      dum2 = index{ncfg,2};
      dum2(end+1:30) = ' ';
    end
  end
end

index = sortrows(index(:,[2 3 1]));
index = index(:, [3 1 2]);
count = 0;
for i=2:size(index,1)
  prevfun = index{i-1,1};
  prevcfg = index{i-1,2};
  prevcmt = index{i-1,3};
  thisfun = index{i,1};
  thiscfg = index{i,2};
  thiscmt = index{i,3};

  if strcmp(thiscfg,prevcfg) && strcmp(thiscmt,prevcmt)
    count = count + 1;
    thisfun = [prevfun ', ' thisfun];
    prevfun = '';
    prevcfg = '';
    prevcmt = '';
    index{i  ,1} = thisfun;
    index{i-1,1} = prevfun;
    index{i-1,2} = prevcfg;
    index{i-1,3} = prevcmt;
  end
end
fprintf('merged %d cfg options\n', count);

fid = fopen(filename, 'wb');
currletter = char(96);
fprintf(fid, '====== Index of configuration options ======\n\n');
fprintf(fid, 'A detailed description of each function is available in the [[:reference|reference documentation]].\n\n');
fprintf(fid, 'This index to the reference documentation is automatically generated from the Matlab code every day. Therefore you should not edit this pages manually, since your changes would be overwritten automatically. If you want to suggest corrections to the documentation, please send them by email to the mailing list or to one of the main developers (see [[:contact]]).\n\n');
for i=1:size(index,1)
  fprintf('%s -- %s -- %s\n', index{i,2}, index{i,3}, index{i,1});
  if isempty(index{i,1})
    continue;
  elseif length(index{i,2})<5
    continue;
  end
  thisletter = index{i,2}(5);
  while currletter<thisletter
    currletter = currletter + 1;
    fprintf(fid, '===== %s =====\n\n', upper(char(currletter)));
  end
  fprintf(fid, '** %s ** // %s //\\\\\n', index{i,2}, index{i,1});
  fprintf(fid, '%s\n\n', index{i,3});
end
fclose(fid);

return

