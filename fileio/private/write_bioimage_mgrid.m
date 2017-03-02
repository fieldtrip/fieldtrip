function write_bioimage_mgrid(filename, elec)

% --------------------------------------------------------
% WRITE_BIOIMAGE_MGRID writes BioImage Suite .mgrid files from a FieldTrip
% elec datatype structure
%
% Use as:
%   write_bioimage_mgrid(filename, elec)
%   where filename has an .mgrid file extension and elec has both a label
%   and an elecpos field
%
% Copyright (C) 2016, Arjen Stolk & Sandon Griffin
% --------------------------------------------------------


% extract info from label field
for e = 1:numel(elec.label) % electrode loop
  ElecStrs{e,1} = regexprep(elec.label{e}, '\d+(?:_(?=\d))?', ''); % without electrode numbers
  ElecNrs(e) = str2double(regexp(elec.label{e},'-?\d+\.?\d*|-?\d*\.?\d+', 'match')); % without electrode strings
end
GridDescript = unique(ElecStrs);
ngrids = numel(GridDescript);
for g = 1:ngrids % grid loop
  Grid2Elec{g} = match_str(ElecStrs, GridDescript{g}); % assign electrodes to grids
end

% open and write ascii-file line by line
fid = fopen(filename, 'wt'); % open ascii-file

% file header
fprintf(fid, '#vtkpxElectrodeMultiGridSource File\n');
fprintf(fid, '#Description\n');
fprintf(fid, 'patient\n');
fprintf(fid, '#Comment\n');
fprintf(fid, 'no additional comment\n');
fprintf(fid, '#Number of Grids\n');
fprintf(fid, [' ' num2str(ngrids) '\n']);

for g = 1:ngrids % grid loop
  
  % grid info
  fprintf(fid, '#- - - - - - - - - - - - - - - - - - -\n');
  fprintf(fid, ['# Electrode Grid ' num2str(ngrids-1) '\n']); % mgrid count starts at 0
  fprintf(fid, '- - - - - - - - - - - - - - - - - - -\n');
  fprintf(fid, '#vtkpxElectrodeGridSource File v2\n');
  fprintf(fid, '#Description\n');
  fprintf(fid, [GridDescript{g} '\n']);
  fprintf(fid, '#Dimensions\n');
  
  % determine grid dimensions
  if isequal(numel(Grid2Elec{g}), 256)
    GridDim(1) = 16; GridDim(2) = 16;
  elseif isequal(numel(Grid2Elec{g}), 64)
    GridDim(1) = 8; GridDim(2) = 8;
  elseif isequal(numel(Grid2Elec{g}), 32)
    e4 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(4)]),:);
    e5 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(5)]),:);
    d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
    if d4to5 < 15 % within 15 mm
      GridDim(1) = 4; GridDim(2) = 8;
    else
      GridDim(1) = 8; GridDim(2) = 4;
    end
  elseif isequal(numel(Grid2Elec{g}), 48)
    e6 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(6)]),:);
    e7 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(7)]),:);
    d6to7 = sqrt(sum((e6-e7).^2)); % distance of elec 6 to 7
    if d6to7 < 15
      GridDim(1) = 6; GridDim(2) = 8;
    else
      GridDim(1) = 8; GridDim(2) = 6;
    end
  elseif isequal(numel(Grid2Elec{g}), 20)
    e4 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(4)]),:);
    e5 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(5)]),:);
    d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
    if d4to5 < 15
      GridDim(1) = 4; GridDim(2) = 5;
    else
      GridDim(1) = 5; GridDim(2) = 4;
    end
  elseif isequal(numel(Grid2Elec{g}), 16)
    e4 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(4)]),:);
    e5 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(5)]),:);
    e9 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(9)]),:);
    d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
    d4to9 = sqrt(sum((e4-e9).^2)); % distance of elec 4 to 9
    if d4to5 > 15 && d4to9 > 15
      GridDim(1) = 4; GridDim(2) = 4;
    elseif d4to5 < 15 && d4to9 > 15
      GridDim(1) = 2; GridDim(2) = 8;
    else
      GridDim(1) = 1; GridDim(2) = 16;
    end
  elseif isequal(numel(Grid2Elec{g}), 12)
    e6 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(6)]),:);
    e7 = elec.elecpos(match_str(elec.label, [GridDescript{g} num2str(7)]),:);
    d6to7 = sqrt(sum((e6-e7).^2)); % distance of elec 6 to 7
    if d6to7 > 15
      GridDim(1) = 2; GridDim(2) = 6;
    else
      GridDim(1) = 1; GridDim(2) = 12;
    end
  else
    GridDim(1) = 1; GridDim(2) = numel(Grid2Elec{g}); % depth or single row strip (not tested)
  end
  fprintf(fid, [' ' num2str(GridDim(1)) ' ' num2str(GridDim(2)) '\n']);
  fprintf(fid, '#Electrode Spacing\n');
  fprintf(fid, ' 10.0000 10.0000\n');
  fprintf(fid, '#Electrode Type\n');
  fprintf(fid, '0\n');
  fprintf(fid, '#Radius\n');
  fprintf(fid, '2.000000\n');
  fprintf(fid, '#Thickeness\n');
  fprintf(fid, '0.050000\n');
  fprintf(fid, '#Color\n');
  fprintf(fid, [num2str(rand) ' ' num2str(rand) ' ' num2str(rand) '\n']);
  
  for e = 1:numel(Grid2Elec{g}) % elec loop
    
    % electrode info
    fprintf(fid, '#- - - - - - - - - - - - - - - - - - -\n');
    ElecNr(1) = mod(e-1,GridDim(1))+1; % 1, 2, .. GridDim, 1, 2, ..
    ElecNr(2) = ceil(e/GridDim(1)); % 1, 2, .. Inf
    fprintf(fid, ['# Electrode ' num2str(ElecNr(1)-1) ' ' num2str(ElecNr(2)-1) '\n']); % mgrid count starts at 0
    fprintf(fid, '- - - - - - - - - - - - - - - - - - -\n');
    fprintf(fid, '#vtkpxElectrodeSource2 File\n');
    fprintf(fid, '#Position\n');
    fprintf(fid, [' ' num2str(elec.elecpos(Grid2Elec{g}(e),1)) ' ' ...
      num2str(elec.elecpos(Grid2Elec{g}(e),2)) ' ' num2str(elec.elecpos(Grid2Elec{g}(e),3)) '\n']);
    fprintf(fid, '#Normal\n');
    fprintf(fid, ' 1.0000 0.0000 0.0000\n');
    fprintf(fid, '#Motor Function\n');
    fprintf(fid, '0\n');
    fprintf(fid, '#Sensory Function\n');
    fprintf(fid, '0\n');
    fprintf(fid, '#Visual Function\n');
    fprintf(fid, '0\n');
    fprintf(fid, '#Language Function\n');
    fprintf(fid, '0\n');
    fprintf(fid, '#Auditory Function\n');
    fprintf(fid, '0\n');
    fprintf(fid, '#User1 Function\n');
    fprintf(fid, '0\n');
    fprintf(fid, '#User2 Function\n');
    fprintf(fid, '0\n');
    fprintf(fid, '#Seizure Onset\n');
    fprintf(fid, '0\n');
    fprintf(fid, '#Spikes Present\n');
    fprintf(fid, '0\n');
    fprintf(fid, '#Electrode Present\n');
    fprintf(fid, '1\n');
    fprintf(fid, '#Electrode Type\n');
    fprintf(fid, '0\n');
    fprintf(fid, '#Radius\n');
    fprintf(fid, '2.000000\n');
    fprintf(fid, '#Thickeness\n');
    fprintf(fid, '0.050000\n');
    fprintf(fid, '#Values\n');
    fprintf(fid, '1\n');
    fprintf(fid, '0.000000\n');
    
  end % end of elec loop
end % end of grid loop
fclose(fid);
