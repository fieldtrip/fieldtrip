function X2=cpd_deform(X, sampling, power, sigma, lambda);

if ~exist('sigma','var') || isempty(sigma), sigma = 1; end;
if ~exist('lambda','var') || isempty(lambda), lambda = 1; end;

[n, d]=size(X);

minx=min(X);
maxx=max(X);

if d==2
    [x, y]=meshgrid(minx(1):sampling*(maxx(1)-minx(1)):maxx(1), minx(2):sampling*(maxx(2)-minx(2)):maxx(2));
    Y=[x(:) y(:)];
else
    [x, y, z]=meshgrid(minx(1):sampling*(maxx(1)-minx(1)):maxx(1), minx(2):sampling*(maxx(2)-minx(2)):maxx(2), minx(3):sampling*(maxx(3)-minx(3)):maxx(3));
    Y=[x(:) y(:) z(:)];
end 
[m, d]=size(Y);


Y2=Y+power*randn(size(Y));


U=Y2-Y;
W=(lambda*eye(m)+cpd_G(Y,Y2,sigma)) \ U;
   
X2=X+cpd_G(X,Y,sigma)*W;
