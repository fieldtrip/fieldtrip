function test_bug3475

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_plot_sens

%%
elec = ft_read_sens('standard_1020.elc');

figure
% electrode labels should appear ~1cm above the surface
ft_plot_sens(elec, 'label', 'on')
xlabel('x');
ylabel('y');
zlabel('z');
view(3)

%%

elec = [];
elec.label = {'elec1', 'elec2', 'elec3', 'elec4', 'elec5'}';
elec.elecpos(:,1) = [10 20 30 40 50];
elec.elecpos(:,2) = [1 1 1 1 1];
elec.elecpos(:,3) = [1 1 1 1 1];
elec.chanpos = elec.elecpos;
elec.unit = 'mm';

figure
ft_plot_sens(elec, 'label', 'on')
xlabel('x');
ylabel('y');
zlabel('z');
view(3)


figure
ft_plot_sens(elec, 'label', 'on', 'elecshape', 'point', 'elecsize', 5)
ft_plot_sens(elec, 'label', 'on', 'elecshape', 'circle', 'elecsize', 5)
ft_plot_sens(elec, 'label', 'on', 'elecshape', 'square', 'elecsize', 5)
ft_plot_sens(elec, 'label', 'on', 'elecshape', 'sphere', 'elecsize', 5)
axis on
xlabel('x');
ylabel('y');
zlabel('z');
view(3)

%%

elec = [];
elec.label = {'elec1'}';
elec.elecpos(1,:) = [0 0 0];
elec.chanpos = elec.elecpos;
elec.unit = 'mm';

figure
ft_plot_sens(elec, 'label', 'on')
xlabel('x');
ylabel('y');
zlabel('z');
view(3)
