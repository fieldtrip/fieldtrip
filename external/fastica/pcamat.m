function [E, D] = pcamat(vectors, firstEig, lastEig, s_interactive, ...
    s_verbose);
%PCAMAT - Calculates the pca for data
%
% [E, D] = pcamat(vectors, firstEig, lastEig, ... 
%                 interactive, verbose);
%
% Calculates the PCA matrices for given data (row) vectors. Returns
% the eigenvector (E) and diagonal eigenvalue (D) matrices containing the
% selected subspaces. Dimensionality reduction is controlled with
% the parameters 'firstEig' and 'lastEig' - but it can also be done
% interactively by setting parameter 'interactive' to 'on' or 'gui'.
%
% ARGUMENTS
%
% vectors       Data in row vectors.
% firstEig      Index of the largest eigenvalue to keep.
%               Default is 1.
% lastEig       Index of the smallest eigenvalue to keep.
%               Default is equal to dimension of vectors.
% interactive   Specify eigenvalues to keep interactively. Note that if
%               you set 'interactive' to 'on' or 'gui' then the values
%               for 'firstEig' and 'lastEig' will be ignored, but they
%               still have to be entered. If the value is 'gui' then the
%               same graphical user interface as in FASTICAG will be
%               used. Default is 'off'.
% verbose       Default is 'on'.
%
%
% EXAMPLE
%       [E, D] = pcamat(vectors);
%
% Note 
%       The eigenvalues and eigenvectors returned by PCAMAT are not sorted.
%
% This function is needed by FASTICA and FASTICAG

% For historical reasons this version does not sort the eigenvalues or
% the eigen vectors in any ways. Therefore neither does the FASTICA or
% FASTICAG. Generally it seams that the components returned from
% whitening is almost in reversed order. (That means, they usually are,
% but sometime they are not - depends on the EIG-command of matlab.)

% @(#)$Id: pcamat.m,v 1.5 2003/12/15 18:24:32 jarmo Exp $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values:
if nargin < 5, s_verbose = 'on'; end
if nargin < 4, s_interactive = 'off'; end
if nargin < 3, lastEig = size(vectors, 1); end
if nargin < 2, firstEig = 1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the optional parameters;
switch lower(s_verbose)
 case 'on'
  b_verbose = 1;
 case 'off'
  b_verbose = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

switch lower(s_interactive)
 case 'on'
  b_interactive = 1;
 case 'off'
  b_interactive = 0;
 case 'gui'
  b_interactive = 2;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''interactive''\n', ...
		s_interactive));
end

oldDimension = size (vectors, 1);
if ~(b_interactive)
  if lastEig < 1 | lastEig > oldDimension
    error(sprintf('Illegal value [ %d ] for parameter: ''lastEig''\n', lastEig));
  end
  if firstEig < 1 | firstEig > lastEig
    error(sprintf('Illegal value [ %d ] for parameter: ''firstEig''\n', firstEig));
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate PCA

% Calculate the covariance matrix.
if b_verbose, fprintf ('Calculating covariance...\n'); end
covarianceMatrix = cov(vectors', 1);

% Calculate the eigenvalues and eigenvectors of covariance
% matrix.
[E, D] = eig (covarianceMatrix);

% The rank is determined from the eigenvalues - and not directly by
% using the function rank - because function rank uses svd, which
% in some cases gives a higher dimensionality than what can be used
% with eig later on (eig then gives negative eigenvalues).
rankTolerance = 1e-7;
maxLastEig = sum (diag (D) > rankTolerance);
if maxLastEig == 0,
  fprintf (['Eigenvalues of the covariance matrix are' ...
	    ' all smaller than tolerance [ %g ].\n' ...
	    'Please make sure that your data matrix contains' ...
	    ' nonzero values.\nIf the values are very small,' ...
	    ' try rescaling the data matrix.\n'], rankTolerance);
  error ('Unable to continue, aborting.');
end

% Sort the eigenvalues - decending.
eigenvalues = flipud(sort(diag(D)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interactive part - command-line
if b_interactive == 1

  % Show the eigenvalues to the user
  hndl_win=figure;
  bar(eigenvalues);
  title('Eigenvalues');

  % ask the range from the user...
  % ... and keep on asking until the range is valid :-)
  areValuesOK=0;
  while areValuesOK == 0
    firstEig = input('The index of the largest eigenvalue to keep? (1) ');
    lastEig = input(['The index of the smallest eigenvalue to keep? (' ...
                    int2str(oldDimension) ') ']);
    % Check the new values...
    % if they are empty then use default values
    if isempty(firstEig), firstEig = 1;end
    if isempty(lastEig), lastEig = oldDimension;end
    % Check that the entered values are within the range
    areValuesOK = 1;
    if lastEig < 1 | lastEig > oldDimension
      fprintf('Illegal number for the last eigenvalue.\n');
      areValuesOK = 0;
    end
    if firstEig < 1 | firstEig > lastEig
      fprintf('Illegal number for the first eigenvalue.\n');
      areValuesOK = 0;
    end
  end
  % close the window
  close(hndl_win);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interactive part - GUI
if b_interactive == 2

  % Show the eigenvalues to the user
  hndl_win = figure('Color',[0.8 0.8 0.8], ...
    'PaperType','a4letter', ...
    'Units', 'normalized', ...
    'Name', 'FastICA: Reduce dimension', ...
    'NumberTitle','off', ...
    'Tag', 'f_eig');
  h_frame = uicontrol('Parent', hndl_win, ...
    'BackgroundColor',[0.701961 0.701961 0.701961], ...
    'Units', 'normalized', ...
    'Position',[0.13 0.05 0.775 0.17], ...
    'Style','frame', ...
    'Tag','f_frame');

b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'HorizontalAlignment','left', ...
	'Position',[0.142415 0.0949436 0.712077 0.108507], ...
	'String','Give the indices of the largest and smallest eigenvalues of the covariance matrix to be included in the reduced data.', ...
	'Style','text', ...
	'Tag','StaticText1');
e_first = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback',[ ...
          'f=round(str2num(get(gcbo, ''String'')));' ...
          'if (f < 1), f=1; end;' ...
          'l=str2num(get(findobj(''Tag'',''e_last''), ''String''));' ...
          'if (f > l), f=l; end;' ...
          'set(gcbo, ''String'', int2str(f));' ...
          ], ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','right', ...
	'Position',[0.284831 0.0678168 0.12207 0.0542535], ...
	'Style','edit', ...
        'String', '1', ...
	'Tag','e_first');
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'HorizontalAlignment','left', ...
	'Position',[0.142415 0.0678168 0.12207 0.0542535], ...
	'String','Range from', ...
	'Style','text', ...
	'Tag','StaticText2');
e_last = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback',[ ...
          'l=round(str2num(get(gcbo, ''String'')));' ...
          'lmax = get(gcbo, ''UserData'');' ...
          'if (l > lmax), l=lmax; fprintf([''The selected value was too large, or the selected eigenvalues were close to zero\n'']); end;' ...
          'f=str2num(get(findobj(''Tag'',''e_first''), ''String''));' ...
          'if (l < f), l=f; end;' ...
          'set(gcbo, ''String'', int2str(l));' ...
          ], ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','right', ...
	'Position',[0.467936 0.0678168 0.12207 0.0542535], ...
	'Style','edit', ...
        'String', int2str(maxLastEig), ...
        'UserData', maxLastEig, ...
	'Tag','e_last');
% in the first version oldDimension was used instead of 
% maxLastEig, but since the program would automatically
% drop the eigenvalues afte maxLastEig...
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'HorizontalAlignment','left', ...
	'Position',[0.427246 0.0678168 0.0406901 0.0542535], ...
	'String','to', ...
	'Style','text', ...
	'Tag','StaticText3');
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback','uiresume(gcbf)', ...
	'Position',[0.630697 0.0678168 0.12207 0.0542535], ...
	'String','OK', ...
	'Tag','Pushbutton1');
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback',[ ...
          'gui_help(''pcamat'');' ...
          ], ...
	'Position',[0.767008 0.0678168 0.12207 0.0542535], ...
	'String','Help', ...
	'Tag','Pushbutton2');

  h_axes = axes('Position' ,[0.13 0.3 0.775 0.6]);
  set(hndl_win, 'currentaxes',h_axes);
  bar(eigenvalues);
  title('Eigenvalues');

  uiwait(hndl_win);
  firstEig = str2num(get(e_first, 'String'));
  lastEig = str2num(get(e_last, 'String'));

  % close the window
  close(hndl_win);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See if the user has reduced the dimension enought

if lastEig > maxLastEig
  lastEig = maxLastEig;
  if b_verbose
    fprintf('Dimension reduced to %d due to the singularity of covariance matrix\n',...
           lastEig-firstEig+1);
  end
else
  % Reduce the dimensionality of the problem.
  if b_verbose
    if oldDimension == (lastEig - firstEig + 1)
      fprintf ('Dimension not reduced.\n');
    else
      fprintf ('Reducing dimension...\n');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drop the smaller eigenvalues
if lastEig < oldDimension
  lowerLimitValue = (eigenvalues(lastEig) + eigenvalues(lastEig + 1)) / 2;
else
  lowerLimitValue = eigenvalues(oldDimension) - 1;
end

lowerColumns = diag(D) > lowerLimitValue;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drop the larger eigenvalues
if firstEig > 1
  higherLimitValue = (eigenvalues(firstEig - 1) + eigenvalues(firstEig)) / 2;
else
  higherLimitValue = eigenvalues(1) + 1;
end
higherColumns = diag(D) < higherLimitValue;

% Combine the results from above
selectedColumns = lowerColumns & higherColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print some info for the user
if b_verbose
  fprintf ('Selected [ %d ] dimensions.\n', sum (selectedColumns));
end
if sum (selectedColumns) ~= (lastEig - firstEig + 1),
  error ('Selected a wrong number of dimensions.');
end

if b_verbose
  fprintf ('Smallest remaining (non-zero) eigenvalue [ %g ]\n', eigenvalues(lastEig));
  fprintf ('Largest remaining (non-zero) eigenvalue [ %g ]\n', eigenvalues(firstEig));
  fprintf ('Sum of removed eigenvalues [ %g ]\n', sum(diag(D) .* ...
    (~selectedColumns)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select the colums which correspond to the desired range
% of eigenvalues.
E = selcol(E, selectedColumns);
D = selcol(selcol(D, selectedColumns)', selectedColumns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some more information
if b_verbose
  sumAll=sum(eigenvalues);
  sumUsed=sum(diag(D));
  retained = (sumUsed / sumAll) * 100;
  fprintf('[ %g ] %% of (non-zero) eigenvalues retained.\n', retained);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newMatrix = selcol(oldMatrix, maskVector);

% newMatrix = selcol(oldMatrix, maskVector);
%
% Selects the columns of the matrix that marked by one in the given vector.
% The maskVector is a column vector.

% 15.3.1998

if size(maskVector, 1) ~= size(oldMatrix, 2),
  error ('The mask vector and matrix are of uncompatible size.');
end

numTaken = 0;

for i = 1 : size (maskVector, 1),
  if maskVector(i, 1) == 1,
    takingMask(1, numTaken + 1) = i;
    numTaken = numTaken + 1;
  end
end

newMatrix = oldMatrix(:, takingMask);