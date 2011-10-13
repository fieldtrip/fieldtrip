function test_ft_average_sens

% TEST test_ft_average_sens
% TEST ft_average_sens ft_transform_sens ft_transform_geometry

% this script so far only tests whether it runs through for a set of eeg
% or meg sensor arrays

% create a set of dummy sensor arrays
sens1 = [];
sens1.pnt   = randn(10,3);
sens1.label = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'}';

sens2 = [];
sens2.pnt   = randn(10,3);
sens2.label = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'}';

sens(1) = sens1;
sens(2) = sens2;

eeg = ft_average_sens(sens);

% also do this for meg
clear sens;

sens1 = [];
sens1.pnt   = randn(10,3);
sens1.ori   = randn(10,3);
sens1.tra   = eye(10);
sens1.label = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'}';

sens2 = [];
sens2.pnt   = randn(10,3);
sens2.ori   = randn(10,3);
sens2.tra   = eye(10);
sens2.label = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'}';

sens(1) = sens1;
sens(2) = sens2;

meg = ft_average_sens(sens);
