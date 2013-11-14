function P = sammon(D, P, varargin)

%SAMMON Computes Sammon's mapping of a data set.
%
% P = sammon(D, P, [value], [mode], [alpha], [Mdist])
%
%  P = sammon(D,2);            % projection to 2-dim space
%  P = sammon(sMap,3);         % projects the codebook vectors
%  P = sammon(sMap,3,[],[],[],Md) % uses distance matrix Md
%  som_grid(sMap,'Coord',P)    % visualization of map projection
%
%  Input and output arguments ([]'s are optional):
%   D        (matrix) size dlen x dim, data to be projected
%            (struct) data or map struct            
%   P        (scalar) output dimension
%            (matrix) size dlen x odim, initial projection matrix
%   [value]  (scalar) all different modes (the next argument) require 
%                     a value, default = 100
%   [mode]   (string) 'steps' or 'errlimit' or 'errchange' or 'seconds',
%                     see below, default is 'steps'
%   [alpha]  (scalar) iteration step size, default = 0.2
%   [Dist]   (matrix) pairwise distance matrix, size dlen x dlen.
%                     If the distances in the input space should
%                     be calculated otherwise than as euclidian
%                     distances, the distance from each vector
%                     to each other vector can be given here,
%                     size dlen x dlen. For example PDIST
%                     function can be used to calculate the
%                     distances: Dist = squareform(pdist(D,'mahal'));
%
%   P        (matrix) size dlen x odim, the projections
%
% The output dimension must be 2 or higher but (naturally) lower 
% than data set dimension.
%
% The mode argument determines the end condition for iteration. If 
% the mode argument is used, also the value argument has to be 
% specified. Different mode possibilities are:
% 'steps'      the iteration is terminated when it is run <value> 
% 'errlimit'   steps, the iteration is terminated when projection error 
%              is lower than <value>,
% 'errchange'  the iteration is terminated when change between 
%              projection error on two successive iteration rounds
%	       is less than <value> percent of total error, and
% 'seconds'    the iteration is terminated after <value> seconds 
%              of iteration.
%
% See also CCA, PCAPROJ, SOM_GRID.

% Reference: Sammon, J.W. Jr., "A nonlinear mapping for data
%   structure analysis", IEEE Transactions on Computers, vol. C-18,
%   no. 5, 1969, pp. 401-409.

% Contributed to SOM Toolbox vs2, February 2nd, 2000 by Juha Vesanto
% Copyright (c) by Juha Vesanto
% http://www.cis.hut.fi/projects/somtoolbox/

% juuso 040100 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check arguments

error(nargchk(2, 6, nargin));  % check no. of input arguments is correct

% input data
if isstruct(D),
  if isfield(D, 'data'),         D = D.data; % data struct
  elseif isfield(D, 'codebook'), D = D.codebook; % map struct
  else error('Invalid structure');
  end
end
if any(isnan(D(:))), 
  error('Cannot make Sammon''s projection for data with unknown components')
end 

% compute data dimensions
orig_si = size(D); 
dim = orig_si(end); 
noc = prod(orig_si)/dim;
if length(orig_si)>2, D = reshape(D,[noc dim]); end

% output dimension / initial projection matrix
if prod(size(P))==1, 
  odim = P; 
  P = rand(noc,odim)-0.5; 
else 
  si = size(P);
  odim = si(end);
  if prod(si) ~= noc*odim, 
    error('Initial projection matrix size does not match data size');
  end
  if length(si)>2, P = reshape(P,[noc odim]); end
  inds = find(isnan(P)); 
  if length(inds), P(inds) = rand(size(inds)); end
end
if odim > dim | odim < 2, 
  error('Output dimension must be within [2, dimension of data]');
end

% determine operating mode
if nargin < 3 | isempty(varargin{1}) | isnan(varargin{1}), value=100;
else value = varargin{1}; 
end
  
if nargin < 4 | isempty(varargin{2}) | isnan(varargin{2}), mode='steps';
else mode = varargin{2}; 
end  
switch mode,
case 'steps',     runlen = value;
case 'errlimit',  errlimit = value;
case 'errchange', errchange = value; e_prev = 0;
case 'seconds',   endtime = value;
otherwise, error(['Illegal mode: ' mode]);
end

% iteration step size
if nargin > 4, alpha = varargin{3}; else alpha = NaN; end
if isempty(alpha) | isnan(alpha), alpha = 0.2; end

% mutual distances
if nargin > 5, Mdist = varargin{4}; else Mdist = []; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialization

% these are used quite frequently
noc_x_1  = ones(noc, 1); 
odim_x_1 = ones(odim,1); 

% compute mutual distances between vectors
if isempty(Mdist) | all(isnan(Mdist(:))),  
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% action

if strcmp(mode, 'seconds'), tic; end;
fprintf(2, 'iterating                    \r');

% sammon iteration   

x  = P ;
xu = zeros(noc, odim);
xd = zeros(noc, odim);
dq = zeros(noc, 1);
dr = zeros(noc, 1);

i = 0;
ready = 0;
while ~ready
  for j = 1:noc,
    xd      = -x + x(j*noc_x_1,:);
    xd2     = xd.^2;
    dpj     = sqrt(sum(xd2'))';
    dq      = Mdist(:,j) - dpj;
    dr      = Mdist(:,j) .* dpj;
    ind     = find(dr ~= 0);
    term    = dq(ind) ./ dr(ind);
    e1      = sum(xd(ind,:) .* term(:,odim_x_1));
    term2   = ((1.0 + dq(ind) ./ dpj(ind)) ./ dpj(ind)) ./ dr(ind);
    e2      = sum(term) - sum(xd2(ind,:) .* term2(:,odim_x_1));
    xu(j,:) = x(j,:) + alpha * e1 ./ abs(e2);
  end

  % move the center of mass to the center 

  c = sum(xu) / noc;
  x = xu - c(noc_x_1, :);

  i = i + 1;

  % compute mapping error
  % doing this adds about 25% to computing time  
  if 0,
    e = 0; tot = 0;
    for j = 2:noc, 
      d   = Mdist(1:(j - 1), j);
      tot = tot + sum(d);
      ind = find(d ~= 0);
      xd  = -x(1:(j - 1), :) + x(j * ones(j - 1, 1), :);
      ee  = d - sqrt(sum(xd'.^2))';
      e   = e + sum(ee(ind).^2 ./ d(ind));
    end
    e = e/tot; 
    fprintf(2, '\r%d iterations, error %f', i, e);
  else
    fprintf(2, '\r%d iterations', i);
  end
  
  % determine is the iteration ready
  
  switch mode
    case 'steps', 
      if i == runlen, ready = 1; end;
    case 'errlimit',
      if e < errlimit, ready = 1; end;
    case 'errchange',
      if i > 1
	change = 100 * abs(e - e_prev) / e_prev;
	if change < errchange, ready = 1; end;
	fprintf(2, ', change of error %f %%    ', change);
      end
      e_prev = e;
    case 'seconds'
      if toc > endtime, ready = 1; end;
      fprintf(2, ', elapsed time %f seconds  ', toc);
  end
  fprintf(2, '        ');
  
  % If you want to see the Sammon's projection plotted (in 2-D and 3-D case),
  % execute the code below; it is not in use by default to speed up 
  % computation.  
  if 0, 
    clf
    if odim == 1,     plot(x(:,1), noc_x_1, 'o');
    elseif odim == 2, plot(x(:,1), x(:,2), 'o');
    else              plot3(x(:,1), x(:,2), x(:,3), 'o')
    end
    drawnow
  end
end

fprintf(2, '\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up

% reshape
orig_si(end) = odim; 
P = reshape(x, orig_si);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%