function test_bug1600

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_checkdata
% TEST ft_datatype_source

% The problem: ft_checkdata(volume, 'datatype', 'source') does not seem to
% convert the inside back to vectorial representation

% verify the problem
volume         = [];
volume.dim     = [3 4 5];
volume.anatomy = zeros(3,4,5);
volume.anatomy(2,2:3,2:4) = rand(1,2,3);
volume.inside  = volume.anatomy>0;
volume.transform = eye(4);

source = ft_checkdata(volume, 'datatype', 'source');

assert(all(size(source.inside)==[60 1]));
