function failed_ft_datatype_sens

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_ft_datatype_sens
% TEST ft_datatype_sens ft_read_sens ft_read_header ctf2grad bti2grad itab2grad fif2grad yokogawa2grad


%% define the test data
path = dccnpath('/home/common/matlab/fieldtrip/data/test/original');

filename = {
  fullfile('meg', 'ctf151', 'Subject01.ds')
  % 'ctf64/Wat123r1raw.ds' % this one fails because the ctf64 dataset requires special headerformat/dataformat flags
  fullfile('meg', 'ctf275', 'A0132_Aud-Obj-Recognition_20051115_02.ds')
  % 'neuromag122/nmStim150.fif' % this one is not there anymore or got renamed
  fullfile('meg', 'neuromag122', 'jg_single_01raw.fif') % but this one is there, now
  fullfile('meg', 'neuromag306', 'raw.fif')
  fullfile('meg', 'itab153', 'srgcst85_0105.raw')
  fullfile('meg', 'itab153', 'srgcst85_0105.raw.mhd')
  fullfile('meg', 'yokogawa160', 'Continuous1.con')
  fullfile('meg', 'itab28', 'gibb0101.raw')
  fullfile('meg', 'itab28', 'gibb0101.raw.mhd')
  fullfile('meg', 'bti148', 'c,rfhp0.1Hz')
  fullfile('meg', 'bti248', 'e,rfDC')
  fullfile('meg', 'bti248', 'e,rfDC,F,a')
  fullfile('meg', 'itab28_old', 'gibb0101.raw')
  fullfile('meg', 'itab28_old', 'gibb0101.raw.mhd')
  fullfile('meg', 'yokogawa64', '2011_01_28_0354_ME053_AEF.con')
  fullfile('meg', 'yokogawa440', 'S1_MEG_Epoch.raw')
  fullfile('electrodes', 'asa', 'andrew_1020.elc')
  fullfile('electrodes', 'asa', 'standard_1005.elc')
  fullfile('electrodes', 'asa', 'standard_1020.elc')
  fullfile('electrodes', 'asa', 'standard_alphabetic.elc')
  fullfile('electrodes', 'asa', 'standard_postfixed.elc')
  fullfile('electrodes', 'asa', 'standard_prefixed.elc')
  fullfile('electrodes', 'asa', 'standard_primed.elc')
  };


%% run the test

for i=1:length(filename)
  
  % this is for running it on a mentat node
  dataset = fullfile(path, filename{i});
  
  sens1 = ft_read_sens(dataset);
  
  if ft_senstype(sens1, 'eeg')
    
    assert(isfield(sens1, 'elecpos'), 'elecpos not present');
    assert(isfield(sens1, 'chanpos'), 'chanpos not present');
    assert(isfield(sens1, 'unit'), 'unit not present');
    % The following fields are optional according to the documentation in ft_datatype_sens
    % assert(isfield(sens1, 'type'), 'type not present');
    % assert(isfield(sens1, 'chantype'), 'chantype not present');
    % assert(isfield(sens1, 'chanunit'), 'chanunit not present');
    
    % make a degraded version of the sensor array
    sens2 = [];
    sens2.pnt = sens1.elecpos;
    sens2.label = sens1.label;
    
  else
    
    assert(isfield(sens1, 'coilpos'), 'coilpos not present');
    assert(isfield(sens1, 'coilori'), 'coilori not present');
    assert(isfield(sens1, 'chanpos'), 'chanpos not present');
    assert(isfield(sens1, 'chanori'), 'chanori not present');
    assert(isfield(sens1, 'unit'), 'unit not present');
    assert(isfield(sens1, 'type'), 'type not present'); % in principle it is optional, but it is expected to be there for MEG
    assert(isfield(sens1, 'chantype'), 'chantype not present'); % in principle it is optional, but it is expected to be there for MEG
    assert(isfield(sens1, 'chanunit'), 'chanunit not present'); % in principle it is optional, but it is expected to be there for MEG
    
    % make a degraded version of the sensor array
    sens2 = [];
    sens2.pnt = sens1.coilpos;
    sens2.ori = sens1.coilori;
    sens2.tra = sens1.tra;
    sens2.label = sens1.label;
    if isfield(sens1, 'balance')
      sens2.balance = sens1.balance;
    end
    
  end % if eeg or meg
  
  % reconstruct the representation with chanpos and coilpos, this should also add the chantype and chanunit fields
  sens2 = ft_datatype_sens(sens2);
  
  disp(dataset);
  disp('');
  disp(sens1)
  disp(sens2)
  
  % they should be identical
  if isfield(sens1, 'unit')
    sens2 = ft_convert_units(sens2, sens1.unit);
  end
  
  if isfield(sens1, 'type') && strcmp(sens1.type, 'yokogawa440')
    % these are known to be different in the two versions and it is not clear which one is correct
    warning('removing chanpos and chanori for yokogawa440, these are not consistent');
    sens1 = rmfield(sens1, {'chanpos', 'chanori'});
    sens2 = rmfield(sens2, {'chanpos', 'chanori'});
  end
  
  % remove coordsys field as these were not yet present in reference files
  if isfield(sens1, 'coordsys')
    sens1 = rmfield(sens1, 'coordsys');
  end
  
  assert(isequal(sens1, sens2));
end

