% function test_bug2085

% TEST test_bug2085
% TEST ft_senstype ft_senslabel

%% create a volume conductor
vol = [];
vol.r = 10;
vol.o = [0 0 0];
vol.unit = 'cm';
vol = ft_datatype_headmodel(vol);

%% create a set of sensors
[pnt, tri] = icosahedron162;
pnt = pnt .* 10; % convert to cm
sel = find(pnt(:,3)>0);
grad.pnt = pnt(sel,:) .* 1.2;
grad.ori = pnt(sel,:);
grad.tra = eye(length(sel));
for i=1:length(sel)
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
    grad.ori(i,:) = grad.ori(i,:) ./ norm(grad.ori(i,:));
    grad.label{i} = sprintf('magnetometer%d', i);
=======
  grad.ori(i,:) = grad.ori(i,:) ./ norm(grad.ori(i,:));
  grad.label{i} = sprintf('magnetometer%d', i);
>>>>>>> enhancement - created test script for timing the leadfield computation
=======
    grad.ori(i,:) = grad.ori(i,:) ./ norm(grad.ori(i,:));
    grad.label{i} = sprintf('magnetometer%d', i);
>>>>>>> enhancement - changed the persistent variable handling in ft_senslabel, this should be much faster. Added it to the test script. See bug 2085
=======
  grad.ori(i,:) = grad.ori(i,:) ./ norm(grad.ori(i,:));
  grad.label{i} = sprintf('magnetometer%d', i);
>>>>>>> enhancement - created test script for timing the leadfield computation
=======
    grad.ori(i,:) = grad.ori(i,:) ./ norm(grad.ori(i,:));
    grad.label{i} = sprintf('magnetometer%d', i);
>>>>>>> enhancement - changed the persistent variable handling in ft_senslabel, this should be much faster. Added it to the test script. See bug 2085
end
grad.unit = 'cm';

grad1 = ft_datatype_sens(grad);
grad2 = ft_datatype_sens(grad);

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> enhancement - include explicit chantype in the test script, see #2085
=======
>>>>>>> enhancement - merged the git and svn change to the test script
=======
>>>>>>> enhancement - include explicit chantype in the test script, see #2085
% this makes ft_senstype much faster
grad2.type = 'magnetometer'; 
for i=1:length(sel)
  grad2.chantype{i} = 'mag';
  grad2.chanunit{i} = 'T';
end
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
grad2.type = 'magnetometer'; % this makes ft_senstype much faster
>>>>>>> enhancement - created test script for timing the leadfield computation
=======
>>>>>>> enhancement - include explicit chantype in the test script, see #2085
=======
>>>>>>> enhancement - merged the git and svn change to the test script
=======
grad2.type = 'magnetometer'; % this makes ft_senstype much faster
>>>>>>> enhancement - created test script for timing the leadfield computation
=======
>>>>>>> enhancement - include explicit chantype in the test script, see #2085

%% determine the time to compute some leadfields

cfg = [];
cfg.vol = vol;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
cfg.grid.resolution = 4;
cfg.channel = 'all';

tic
=======
cfg.resolution = 4;
=======
cfg.grid.resolution = 4;
>>>>>>> enhancement - include explicit chantype in the test script, see #2085
=======
cfg.grid.resolution = 4;
>>>>>>> enhancement - merged the git and svn change to the test script
cfg.channel = 'all';

<<<<<<< HEAD
tic 
>>>>>>> enhancement - created test script for timing the leadfield computation
=======
tic
>>>>>>> enhancement - changed the persistent variable handling in ft_senslabel, this should be much faster. Added it to the test script. See bug 2085
=======
cfg.resolution = 4;
=======
cfg.grid.resolution = 4;
>>>>>>> enhancement - include explicit chantype in the test script, see #2085
cfg.channel = 'all';

<<<<<<< HEAD
tic 
>>>>>>> enhancement - created test script for timing the leadfield computation
=======
tic
>>>>>>> enhancement - changed the persistent variable handling in ft_senslabel, this should be much faster. Added it to the test script. See bug 2085
cfg.grad = grad1;
grid1 = ft_prepare_leadfield(cfg);
time1 = toc;

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
tic
=======
tic 
>>>>>>> enhancement - created test script for timing the leadfield computation
=======
tic
>>>>>>> enhancement - changed the persistent variable handling in ft_senslabel, this should be much faster. Added it to the test script. See bug 2085
=======
tic 
>>>>>>> enhancement - created test script for timing the leadfield computation
=======
tic
>>>>>>> enhancement - changed the persistent variable handling in ft_senslabel, this should be much faster. Added it to the test script. See bug 2085
cfg.grad = grad2;
grid2 = ft_prepare_leadfield(cfg);
time2 = toc;

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
%% compare the time that it took for the two grad versions, the 10 seconds and fraction 2 are rather arbitrary

if time1>10
  error('the leadfield computation takes too long');
end

if (time1/time2)>2
    error('the leadfield computation without grad.type takes too long (%d seconds compared to %d seconds)', time1, time2);
end


%% the following section pertains to a change (improvement) I made to the ft_senslabel function which is meant to speed up subsequent calls using a large set of persistent variables

type = {
    'btiref'
    'bti148'
    'bti148_planar'
    'bti248'
    'bti248_planar'
    'ctfref'
    'ctfheadloc'
    'ctf64'
    'ctf151'
    'ctf151_planar'
    'ctf275'
    'ctf275_planar'
    'neuromag122'
    'neuromag122alt'
    'neuromag306'
    'neuromag306alt'
    'eeg1020'
    'eeg1010'
    'eeg1005'
    'ext1020'
    'biosemi64'
    'biosemi128'
    'biosemi256'
    'egi32'
    'egi64'
    'egi128'
    'egi256'
    'itab28'
    'itab153'
    'itab153_planar'
    'yokogawa9'
    'yokogawa64'
    'yokogawa64_planar'
    'yokogawa160'
    'yokogawa160_planar'
    'yokogawa440'
    'yokogawa440_planar'
    'electrode'
    };

for i=1:length(type)
    clear ft_senslabel
    disp(type{i});
    label1 = ft_senslabel(type{i});
    label2 = ft_senslabel(type{i}); % this is from persistent cache
    assert(identical(label1, label2), sprintf('problem in senslabel for %s', type{i}));
end




=======
%% compare the time that it took, this fraction 1.5 is rather arbitrary
=======
%% compare the time that it took for the two grad versions, the 10 seconds and fraction 2 are rather arbitrary
>>>>>>> enhancement - do not deal with backward compatible grad in ft_compute_leadfield, leave this to ft_prepare_vol_sens. Updated test script for #2085.

if time1>10
  error('the leadfield computation takes too long');
end

if (time1/time2)>2
    error('the leadfield computation without grad.type takes too long (%d seconds compared to %d seconds)', time1, time2);
end

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> enhancement - created test script for timing the leadfield computation
=======
%%
=======
=======

>>>>>>> enhancement - do not deal with backward compatible grad in ft_compute_leadfield, leave this to ft_prepare_vol_sens. Updated test script for #2085.
%% the following section pertains to a change (improvement) I made to the ft_senslabel function which is meant to speed up subsequent calls using a large set of persistent variables
>>>>>>> enhancement - include explicit chantype in the test script, see #2085
=======
%% the following section pertains to a change (improvement) I made to the ft_senslabel function which is meant to speed up subsequent calls using a large set of persistent variables
>>>>>>> enhancement - merged the git and svn change to the test script

type = {
    'btiref'
    'bti148'
    'bti148_planar'
    'bti248'
    'bti248_planar'
    'ctfref'
    'ctfheadloc'
    'ctf64'
    'ctf151'
    'ctf151_planar'
    'ctf275'
    'ctf275_planar'
    'neuromag122'
    'neuromag122alt'
    'neuromag306'
    'neuromag306alt'
    'eeg1020'
    'eeg1010'
    'eeg1005'
    'ext1020'
    'biosemi64'
    'biosemi128'
    'biosemi256'
    'egi32'
    'egi64'
    'egi128'
    'egi256'
    'itab28'
    'itab153'
    'itab153_planar'
    'yokogawa9'
    'yokogawa64'
    'yokogawa64_planar'
    'yokogawa160'
    'yokogawa160_planar'
    'yokogawa440'
    'yokogawa440_planar'
    'electrode'
    };

for i=1:length(type)
    clear ft_senslabel
    disp(type{i});
    label1 = ft_senslabel(type{i});
    label2 = ft_senslabel(type{i}); % this is from persistent cache
    assert(identical(label1, label2), sprintf('problem in senslabel for %s', type{i}));
end




>>>>>>> enhancement - changed the persistent variable handling in ft_senslabel, this should be much faster. Added it to the test script. See bug 2085
=======
%% compare the time that it took, this fraction 1.5 is rather arbitrary

if (time1/time2)>1.5
    error('the leadfield computation without grad.type takes too long (%d seconds compared to %d seconds)', time1, time2);
end

<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> enhancement - created test script for timing the leadfield computation
=======
%%
=======
%% the following section pertains to a change (improvement) I made to the ft_senslabel function which is meant to speed up subsequent calls using a large set of persistent variables
>>>>>>> enhancement - include explicit chantype in the test script, see #2085

type = {
    'btiref'
    'bti148'
    'bti148_planar'
    'bti248'
    'bti248_planar'
    'ctfref'
    'ctfheadloc'
    'ctf64'
    'ctf151'
    'ctf151_planar'
    'ctf275'
    'ctf275_planar'
    'neuromag122'
    'neuromag122alt'
    'neuromag306'
    'neuromag306alt'
    'eeg1020'
    'eeg1010'
    'eeg1005'
    'ext1020'
    'biosemi64'
    'biosemi128'
    'biosemi256'
    'egi32'
    'egi64'
    'egi128'
    'egi256'
    'itab28'
    'itab153'
    'itab153_planar'
    'yokogawa9'
    'yokogawa64'
    'yokogawa64_planar'
    'yokogawa160'
    'yokogawa160_planar'
    'yokogawa440'
    'yokogawa440_planar'
    'electrode'
    };

for i=1:length(type)
    clear ft_senslabel
    disp(type{i});
    label1 = ft_senslabel(type{i});
    label2 = ft_senslabel(type{i}); % this is from persistent cache
    assert(identical(label1, label2), sprintf('problem in senslabel for %s', type{i}));
end




>>>>>>> enhancement - changed the persistent variable handling in ft_senslabel, this should be much faster. Added it to the test script. See bug 2085
