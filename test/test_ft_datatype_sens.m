function test_ft_datatype_sens

% TEST test_ft_datatype_sens
% TEST ft_datatype_sens ft_read_sens ft_read_header ctf2grad bti2grad itab2grad fif2grad yokogawa2grad


path1 = '/home/common/matlab/fieldtrip/data/test/original/meg';
path2 = '/Users/robert/Manzana/data/dataformat/testdata';

filename = {
  'ctf151/Subject01.ds'
  % this one fails because the ctf64 dataset requires special headerformat/dataformat flags
  % 'ctf64/Wat123r1raw.ds'
  'ctf275/A0132_Aud-Obj-Recognition_20051115_02.ds'
  'neuromag122/nmStim150.fif'
  'neuromag306/raw.fif'
  'itab153/srgcst85_0105.raw'
  'itab153/srgcst85_0105.raw.mhd'
  'yokogawa160/Continuous1.con'
  'itab28/gibb0101.raw'
  'itab28/gibb0101.raw.mhd'
  'bti148/c,rfhp0.1Hz'
  'bti248/e,rfDC'
  'bti248/e,rfDC,F,a'
  'itab28_old/gibb0101.raw'
  'itab28_old/gibb0101.raw.mhd'
  'yokogawa64/2011_01_28_0354_ME053_AEF.con'
  'yokogawa440/S1_MEG_Epoch.raw'
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
      warning('unable to test with the dataset "%s" because it cannot be found', filename{i});
      continue
    end
  else
    error('unable to test with the dataset "%d" because it cannot be found', filename{i});
  end
  
  sens1 = ft_read_sens(dataset);
  
  assert(isfield(sens1, 'chanpos'));
  assert(isfield(sens1, 'chanori'));
  assert(isfield(sens1, 'chantype'));
  assert(isfield(sens1, 'chanunit'));
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
    assert(isequal(sens1, sens2));
  catch
    warning(lasterr);
  end
end

