function test_ft_datatype_sens

% TEST test_ft_datatype_sens
% TEST ft_datatype_sens ft_read_sens ft_read_header ctf2grad bti2grad itab2grad fif2grad yokogawa2grad

if isunix
  path1 = '/home/common/matlab/fieldtrip/data/test/original/meg';
elseif ispc
  path1 = fullfile('H:', 'common', 'matlab', 'fieldtrip', 'data', 'test', 'original', 'meg');
end
path2 = '/Users/robert/Manzana/data/dataformat/testdata';

filename = {
  fullfile('ctf151', 'Subject01.ds')
  % this one fails because the ctf64 dataset requires special headerformat/dataformat flags
  % 'ctf64/Wat123r1raw.ds'
  fullfile('ctf275', 'A0132_Aud-Obj-Recognition_20051115_02.ds')
  % this one is not there anymore or got renamed
  % fullfile('neuromag122', 'nmStim150.fif')
  % but this one is there, now
  fullfile('neuromag122', 'jg_single_01raw.fif')
  fullfile('neuromag306', 'raw.fif')
  fullfile('itab153', 'srgcst85_0105.raw')
  fullfile('itab153', 'srgcst85_0105.raw.mhd')
  fullfile('yokogawa160', 'Continuous1.con')
  fullfile('itab28', 'gibb0101.raw')
  fullfile('itab28', 'gibb0101.raw.mhd')
  fullfile('bti148', 'c,rfhp0.1Hz')
  fullfile('bti248', 'e,rfDC')
  fullfile('bti248', 'e,rfDC,F,a')
  fullfile('itab28_old', 'gibb0101.raw')
  fullfile('itab28_old', 'gibb0101.raw.mhd')
  fullfile('yokogawa64', '2011_01_28_0354_ME053_AEF.con')
  fullfile('yokogawa440', 'S1_MEG_Epoch.raw')
  };

for i=1:length(filename)
  
  if isdir(path1)
    % this is for running it on a mentat node
    dataset = fullfile(path1, filename{i});
  elseif isdir(path2)
    % the following part allows the code to run on the laptop of Robert
    [p, f, x] = fileparts(filename{i});
    [a, b] = system(sprintf('find %s -name %s', path2, [f x]));
    if length(b)>1
      dataset = deblank(b);
    else
      warning('unable to test with the dataset "%s" because it cannot be found\n', filename{i});
      continue
    end
  else
    error('unable to test with the dataset "%d" because it cannot be found\n', filename{i});
  end
  
  sens1 = ft_read_sens(dataset);
  
  assert(isfield(sens1, 'chanpos'));
  assert(isfield(sens1, 'chanori'));
  assert(isfield(sens1, 'chantype'));
  % assert(isfield(sens1, 'chanunit'));
  % assert(isfield(sens1, 'type'));
  % assert(isfield(sens1, 'unit'));
  
  sens2         = [];
  sens2.pnt     = sens1.coilpos;
  sens2.ori     = sens1.coilori;
  sens2.tra     = sens1.tra;
  sens2.label   = sens1.label;
  if isfield(sens1, 'balance')
    sens2.balance = sens1.balance;
  end
  
  % reconstruct the representation with chanpos and coilpos
  % this should also add the chantype and chanunit fields
  sens2 = ft_datatype_sens(sens2);
  
  disp(dataset);
  disp('');
  disp(sens1)
  disp(sens2)
  
  try
    % they should be identical
    if isfield(sens1, 'unit')
      sens2 = ft_convert_units(sens2, sens1.unit);
    end
    assert(isequal(sens1, sens2));
  catch
    warning(lasterr);
  end
end

