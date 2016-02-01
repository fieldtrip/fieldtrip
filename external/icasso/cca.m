function [P] = cca(D, P, epochs, Mdist, alpha0, lambda0)

%CCA Projects data vectors using Curvilinear Component Analysis.
%
% P = cca(D, P, epochs, [Dist], [alpha0], [lambda0])
%
%  P = cca(D,2,10);           % projects the given data to a plane
%  P = cca(D,pcaproj(D,2),5); % same, but with PCA initialization
%  P = cca(D, 2, 10, Dist);   % same, but the given distance matrix is used
%  
%  Input and output arguments ([]'s are optional):
%   D          (matrix) the data matrix, size dlen x dim
%              (struct) data or map struct            
%   P          (scalar) output dimension
%              (matrix) size dlen x odim, the initial projection
%   epochs     (scalar) training length
%   [Dist]     (matrix) pairwise distance matrix, size dlen x dlen.
%                       If the distances in the input space should
%                       be calculated otherwise than as euclidian
%                       distances, the distance from each vector
%                       to each other vector can be given here,
%                       size dlen x dlen. For example PDIST
%                       function can be used to calculate the
%                       distances: Dist = squareform(pdist(D,'mahal'));
%   [alpha0]   (scalar) initial step size, 0.5 by default
%   [lambda0]  (scalar) initial radius of influence, 3*max(std(D)) by default
%  
%   P          (matrix) size dlen x odim, the projections
%
% Unknown values (NaN's) in the data: projections of vectors with
% unknown components tend to drift towards the center of the
% projection distribution. Projections of totally unknown vectors are
% set to unknown (NaN).
%
% See also SAMMON, PCAPROJ. 

% Reference: Demartines, P., Herault, J., "Curvilinear Component
%   Analysis: a Self-Organizing Neural Network for Nonlinear
%   Mapping of Data Sets", IEEE Transactions on Neural Networks, 
%   vol 8, no 1, 1997, pp. 148-154.

% Contributed to SOM Toolbox 2.0, February 2nd, 2000 by Juha Vesanto
% Copyright (c) by Juha Vesanto
% http://www.cis.hut.fi/projects/somtoolbox/

% juuso 171297 040100
% johan 170305: cca used error as a variable name which caused a
% warning, changed to error_
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check arguments 

error(nargchk(3, 6, nargin)); % check the number of input arguments

% input data
if isstruct(D), 
  if strcmp(D.type,'som_map'), D = D.codebook; else D = D.data; end
end
[noc dim] = size(D);
noc_x_1  = ones(noc, 1); % used frequently
me = zeros(1,dim); st = zeros(1,dim);
for i=1:dim,
  me(i) = mean(D(find(isfinite(D(:,i))),i));
  st(i) = std(D(find(isfinite(D(:,i))),i));
end

% initial projection
if prod(size(P))==1, 
  P = (2*rand(noc,P)-1).*st(noc_x_1,1:P) + me(noc_x_1,1:P); 
else
  % replace unknown projections with known values
  inds = find(isnan(P)); P(inds) = rand(size(inds));
end
[dummy odim] = size(P);
odim_x_1  = ones(odim, 1); % this is used frequently

% training length
train_len = epochs*noc;

% random sample order
rand('state',sum(100*clock));
sample_inds = ceil(noc*rand(train_len,1));

% mutual distances
if nargin<4 | isempty(Mdist) | all(isnan(Mdist(:))),
  fprintf(2, 'computing mutual distances\r');
  dim_x_1 = ones(dim,1);
  for i = 1:noc,
    x = D(i,:); 
    Diff = D - x(noc_x_1,:);
    N = isnan(Diff);
    Diff(find(N)) = 0; 
    Mdist(:,i) = sqrt((Diff.^2)*dim_x_1);
    N = find(sum(N')==dim); %mutual distance unknown
    if ~isempty(N), Mdist(N,i) = NaN; end
  end
else
  % if the distance matrix is output from PDIST function
  if size(Mdist,1)==1, Mdist = squareform(Mdist); end
  if size(Mdist,1)~=noc, 
    error('Mutual distance matrix size and data set size do not match'); 
  end
end

% alpha and lambda
if nargin<5 | isempty(alpha0) | isnan(alpha0), alpha0 = 0.5; end
alpha = potency_curve(alpha0,alpha0/100,train_len);

if nargin<6 | isempty(lambda0) | isnan(lambda0), lambda0 = max(st)*3; end
lambda = potency_curve(lambda0,0.01,train_len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Action

k=0; fprintf(2, 'iterating: %d / %d epochs\r',k,epochs);

for i=1:train_len, 
  
  ind = sample_inds(i);     % sample index
  dx = Mdist(:,ind);        % mutual distances in input space
  known = find(~isnan(dx)); % known distances

  if ~isempty(known),
    % sample vector's projection
    y = P(ind,:);                 

    % distances in output space
    Dy = P(known,:) - y(noc_x_1(known),:); 
    dy = sqrt((Dy.^2)*odim_x_1);           
  
    % relative effect
    dy(find(dy==0)) = 1;        % to get rid of div-by-zero's
    fy = exp(-dy/lambda(i)) .* (dx(known) ./ dy - 1);

    % Note that the function F here is e^(-dy/lambda)) 
    % instead of the bubble function 1(lambda-dy) used in the 
    % paper.
    
    % Note that here a simplification has been made: the derivatives of the
    % F function have been ignored in calculating the gradient of error
    % function w.r.t. to changes in dy.
    
    % update
    P(known,:) = P(known,:) + alpha(i)*fy(:,odim_x_1).*Dy;
  end

  % track
  if rem(i,noc)==0, 
    k=k+1; fprintf(2, 'iterating: %d / %d epochs\r',k,epochs);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear up

% calculate error
error_ = cca_error(P,Mdist,lambda(train_len));
fprintf(2,'%d iterations, error %f          \n', epochs, error_);

% set projections of totally unknown vectors as unknown
unknown = find(sum(isnan(D)')==dim);
P(unknown,:) = NaN;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tips

% to plot the results, use the code below

%subplot(2,1,1), 
%switch(odim), 
%  case 1, plot(P(:,1),ones(dlen,1),'x')
%  case 2, plot(P(:,1),P(:,2),'x'); 
%  otherwise, plot3(P(:,1),P(:,2),P(:,3),'x'); rotate3d on
%end
%subplot(2,1,2), dydxplot(P,Mdist);

% to a project a new point x in the input space to the output space
% do the following:

% Diff = D - x(noc_x_1,:); Diff(find(isnan(Diff))) = 0; 
% dx = sqrt((Diff.^2)*dim_x_1);
% p = project_point(P,x,dx); % this function can be found from below
% tlen = size(p,1);
% plot(P(:,1),P(:,2),'bx',p(tlen,1),p(tlen,2),'ro',p(:,1),p(:,2),'r-')

% similar trick can be made to the other direction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunctions

function vals = potency_curve(v0,vn,l)

  % curve that decreases from v0 to vn with a rate that is 
  % somewhere between linear and 1/t
  vals = v0 * (vn/v0).^([0:(l-1)]/(l-1));


function error_ = cca_error(P,Mdist,lambda)

  [noc odim] = size(P);
  noc_x_1 = ones(noc,1);
  odim_x_1 = ones(odim,1);

  error_ = 0;
  for i=1:noc,
    known = find(~isnan(Mdist(:,i)));
    if ~isempty(known),   
      y = P(i,:);                 
      Dy = P(known,:) - y(noc_x_1(known),:);
      dy = sqrt((Dy.^2)*odim_x_1);
      fy = exp(-dy/lambda);
      error_ = error_ + sum(((Mdist(known,i) - dy).^2).*fy);
    end
  end
  error_ = error_/2;


function [] = dydxplot(P,Mdist)

  [noc odim] = size(P);
  noc_x_1 = ones(noc,1);
  odim_x_1 = ones(odim,1);
  Pdist = zeros(noc,noc);
    
  for i=1:noc,
    y = P(i,:);                 
    Dy = P - y(noc_x_1,:);
    Pdist(:,i) = sqrt((Dy.^2)*odim_x_1);
  end

  Pdist = tril(Pdist,-1); 
  inds = find(Pdist > 0); 
  n = length(inds);
  plot(Pdist(inds),Mdist(inds),'.');
  xlabel('dy'), ylabel('dx')


function p = project_point(P,x,dx)

  [noc odim] = size(P);
  noc_x_1 = ones(noc,1);
  odim_x_1 = ones(odim,1);

  % initial projection
  [dummy,i] = min(dx);
  y = P(i,:)+rand(1,odim)*norm(P(i,:))/20;
 
  % lambda 
  lambda = norm(std(P));

  % termination
  eps = 1e-3; i_max = noc*10;
  
  i=1; p(i,:) = y; 
  ready = 0;
  while ~ready,

    % mutual distances
    Dy = P - y(noc_x_1,:);        % differences in output space
    dy = sqrt((Dy.^2)*odim_x_1);  % distances in output space
    f = exp(-dy/lambda);
  
    fprintf(2,'iteration %d, error %g \r',i,sum(((dx - dy).^2).*f));

    % all the other vectors push the projected one
    fy = f .* (dx ./ dy - 1) / sum(f);
  
    % update    
    step = - sum(fy(:,odim_x_1).*Dy);
    y = y + step;
  
    i=i+1;
    p(i,:) = y;   
    ready = (norm(step)/norm(y) < eps | i > i_max);

  end
  fprintf(2,'\n');
     
  