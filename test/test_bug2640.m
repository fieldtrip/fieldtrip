function test_bug2640

% TEST ft_plot_slice ft_plot_ortho

% WALLTIME 00:10:00
% MEM 1gb

dim = [20, 30, 40];
dat = 0.1*randn(dim);

selx = 5:(dim(1)-5);
sely = 5:(dim(2)-5);
selz = 5:(dim(3)-5);

dat(selx, sely, selz) = dat(selx, sely, selz) + 1;

transform = eye(4);

[f, v] = isosurface(dat>0.5);

mesh.tri = f;
mesh.pnt = v(:,[2 1 3]);

figure
ft_plot_slice(dat, 'transform', transform, 'intersectmesh', mesh);
grid on
xlabel('x')
ylabel('y')
zlabel('z')

figure
ft_plot_ortho(dat, 'transform', transform, 'style', 'intersect', 'intersectmesh', mesh);
grid on
xlabel('x')
ylabel('y')
zlabel('z')

%%%%%%%%%

dim = [3 3 3];
dat = zeros(3,3,3);
dat(2,2,2) = 1;

transform = eye(4);
transform(:,4) = [-2 -2 -2 1]';

mesh.pnt = [
  0 0 0 
  0 0 1
  0 1 0 
  1 0 0 
  0 1 1
  1 0 1
  1 1 0
  1 1 1
  ] - 0.5;
mesh.tri = convhulln(mesh.pnt);

z = linspace(-2.5, 2.5, 20);
figure
for i=1:length(z)
  ft_plot_slice(dat, 'transform', transform, 'intersectmesh', mesh, 'location', [0 0 z(i)]);
end
alpha 0.5
ft_plot_mesh(mesh)

