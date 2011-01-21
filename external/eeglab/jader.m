% jader() - blind separation of real signals using JADE (v1.5, Dec. 1997).
%
% Usage: 
%   >> B = jader(X);
%   >> B = jader(X,m);
%
% Notes:
%  1) If X is an nxT data matrix (n sensors, T samples) then
%     B=jader(X) is a nxn separating matrix such that S=B*X is an nxT
%     matrix of estimated source signals.
%  2) If B=jader(X,m), then B has size mxn so that only m sources are
%     extracted.  This is done by restricting the operation of jader
%     to the m first principal components. 
%  3) Also, the rows of B are ordered such that the columns of pinv(B)
%     are in order of decreasing norm; this has the effect that the
%     `most energetically significant' components appear first in the
%     rows of S=B*X.
%
% Author: Jean-Francois Cardoso (cardoso@sig.enst.fr) 
 
% Quick notes (more at the end of this file)
%
%  o this code is for REAL-valued signals.  An implementation of JADE
%    for both real and complex signals is also available from
%    http://sig.enst.fr/~cardoso/stuff.html
%
%  o This algorithm differs from the first released implementations of
%    JADE in that it has been optimized to deal more efficiently
%    1) with real signals (as opposed to complex)
%    2) with the case when the ICA model does not necessarily hold.
%
%  o There is a practical limit to the number of independent
%    components that can be extracted with this implementation.  Note
%    that the first step of JADE amounts to a PCA with dimensionality
%    reduction from n to m (which defaults to n).  In practice m
%    cannot be `very large' (more than 40, 50, 60... depending on
%    available memory)
%
%  o See more notes, references and revision history at the end of
%    this file and more stuff on the WEB
%    http://sig.enst.fr/~cardoso/stuff.html
%
%  o This code is supposed to do a good job!  Please report any
%    problem to cardoso@sig.enst.fr

% Copyright : Jean-Francois Cardoso.  cardoso@sig.enst.fr

function B =  jadeR(X,m)

verbose	= 1 ;	% Set to 0 for quiet operation


% Finding the number of sources
[n,T]	= size(X);
if nargin==1, m=n ; end; 	% Number of sources defaults to # of sensors
if m>n ,    fprintf('jade -> Do not ask more sources than sensors here!!!\n'), return,end
if verbose, fprintf('jade -> Looking for %d sources\n',m); end ;



% Self-commenting code
%=====================
if verbose, fprintf('jade -> Removing the mean value\n'); end 
X	= X - mean(X')' * ones(1,T);


%%% whitening & projection onto signal subspace
%   ===========================================
if verbose, fprintf('jade -> Whitening the data\n'); end
 [U,D] 		= eig((X*X')/T)	; 
 [puiss,k]	= sort(diag(D))	;
 rangeW		= n-m+1:n			; % indices to the m  most significant directions
 scales		= sqrt(puiss(rangeW))		; % scales
  W  		= diag(1./scales)  * U(1:n,k(rangeW))'	;	% whitener
 iW  		= U(1:n,k(rangeW)) * diag(scales) 	;	% its pseudo-inverse
 X		= W*X;  


%%% Estimation of the cumulant matrices.
%   ====================================
if verbose, fprintf('jade -> Estimating cumulant matrices\n'); end

dimsymm 	= (m*(m+1))/2;	% Dim. of the space of real symm matrices
nbcm 		= dimsymm  ; 	% number of cumulant matrices
CM 		= zeros(m,m*nbcm);  % Storage for cumulant matrices
R 		= eye(m);  	%% 
Qij 		= zeros(m);	% Temp for a cum. matrix
Xim		= zeros(1,m);	% Temp
Xjm		= zeros(1,m);	% Temp
scale		= ones(m,1)/T ; % for convenience



%% I am using a symmetry trick to save storage.  I should write a
%% short note one of these days explaining what is going on here.
%%
Range 		= 1:m ; % will index the columns of CM where to store the cum. mats.
for im = 1:m
	Xim = X(im,:) ;
	Qij = ((scale* (Xim.*Xim)) .* X ) * X' 	- R - 2 * R(:,im)*R(:,im)' ;
	CM(:,Range)	= Qij ; 
	Range 		= Range + m ; 
   for jm = 1:im-1
	Xjm = X(jm,:) ;
	Qij = ((scale * (Xim.*Xjm) ) .*X ) * X' - R(:,im)*R(:,jm)' - R(:,jm)*R(:,im)' ;
	CM(:,Range)	= sqrt(2)*Qij ;  
	Range 		= Range + m ;
   end ;
end;

%%% joint diagonalization of the cumulant matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init
if 1, 	%% Init by diagonalizing a *single* cumulant matrix.  It seems to save
	%% some computation time `sometimes'.  Not clear if initialization is
	%% a good idea since Jacobi rotations are very efficient.

	if verbose, fprintf('jade -> Initialization of the diagonalization\n'); end
	[V,D]	= eig(CM(:,1:m));	% For instance, this one
	for u=1:m:m*nbcm,		% updating accordingly the cumulant set given the init
		CM(:,u:u+m-1) = CM(:,u:u+m-1)*V ; 
	end;
	CM	= V'*CM;

else,	%% The dont-try-to-be-smart init
	V	= eye(m) ; % la rotation initiale
end;

seuil	= 1/sqrt(T)/100; % A statistically significant threshold
encore	= 1;
sweep	= 0;
updates = 0;
g	= zeros(2,nbcm);
gg	= zeros(2,2);
G	= zeros(2,2);
c	= 0 ;
s 	= 0 ;
ton	= 0 ;
toff	= 0 ;
theta	= 0 ;

%% Joint diagonalization proper
if verbose, fprintf('jade -> Contrast optimization by joint diagonalization\n'); end

while encore, encore=0;   

  if verbose, fprintf('jade -> Sweep #%d\n',sweep); end
  sweep=sweep+1;

 for p=1:m-1,
  for q=p+1:m,

 	Ip = p:m:m*nbcm ;
	Iq = q:m:m*nbcm ;

	%%% computation of Givens angle
 	g	= [ CM(p,Ip)-CM(q,Iq) ; CM(p,Iq)+CM(q,Ip) ];
 	gg	= g*g';
 	ton 	= gg(1,1)-gg(2,2); 
 	toff 	= gg(1,2)+gg(2,1);
 	theta	= 0.5*atan2( toff , ton+sqrt(ton*ton+toff*toff) );

 	%%% Givens update
 	if abs(theta) > seuil,	encore = 1 ;
 		updates = updates + 1;
 		c	= cos(theta); 
	 	s	= sin(theta);
 		G	= [ c -s ; s c ] ;

		pair 		= [p;q] ;
		V(:,pair) 	= V(:,pair)*G ;
	 	CM(pair,:)	= G' * CM(pair,:) ;
		CM(:,[Ip Iq]) 	= [ c*CM(:,Ip)+s*CM(:,Iq) -s*CM(:,Ip)+c*CM(:,Iq) ] ;

		%% fprintf('jade -> %3d %3d %12.8f\n',p,q,s);

 	end%%of the if
  end%%of the loop on q
 end%%of the loop on p
end%%of the while loop
if verbose, fprintf('jade -> Total of %d Givens rotations\n',updates); end

%%% A separating matrix
%   ===================
B	= V'*W ;

%%% We permut its rows to get the most energetic components first.
%%% Here the **signals** are normalized to unit variance.  Therefore,
%%% the sort is according to the norm of the columns of A = pinv(B)

if verbose, fprintf('jade -> Sorting the components\n',updates); end
A		= iW*V ;
[vars,keys]	= sort(sum(A.*A)) ;
B		= B(keys,:);
B		= B(m:-1:1,:) ; % Is this smart ?

% Signs are fixed by forcing the first column of B to have
% non-negative entries.
if verbose, fprintf('jade -> Fixing the signs\n',updates); end
b	= B(:,1) ;
signs	= sign(sign(b)+0.1) ; % just a trick to deal with sign=0
B	= diag(signs)*B ;


return ;

% To do.
%   - Implement a cheaper/simpler whitening (is it worth it?)
% 
% Revision history:
%
%-  V1.5, Dec. 24 1997 
%   - The sign of each row of B is determined by letting the first
%     element be positive.  
%
%-  V1.4, Dec. 23 1997 
%   - Minor clean up.
%   - Added a verbose switch
%   - Added the sorting of the rows of B in order to fix in some
%     reasonable way the permutation indetermination.  See note 2)
%     below.
%
%-  V1.3, Nov.  2 1997 
%   - Some clean up.  Released in the public domain.
%
%-  V1.2, Oct.  5 1997 
%   - Changed random picking of the cumulant matrix used for
%     initialization to a deterministic choice.  This is not because
%     of a better rationale but to make the ouput (almost surely)
%     deterministic.
%   - Rewrote the joint diag. to take more advantage of Matlab's
%     tricks.
%   - Created more dummy variables to combat Matlab's loose memory
%     management.
%
%-  V1.1, Oct. 29 1997.
%    Made the estimation of the cumulant matrices more regular. This
%    also corrects a buglet...
%
%-  V1.0, Sept. 9 1997. Created.
%
% Main reference:
% @article{CS-iee-94,
%  title 	= "Blind beamforming for non {G}aussian signals",
%  author       = "Jean-Fran\c{c}ois Cardoso and Antoine Souloumiac",
%  HTML 	= "ftp://sig.enst.fr/pub/jfc/Papers/iee.ps.gz",
%  journal      = "IEE Proceedings-F",
%  month = dec, number = 6, pages = {362-370}, volume = 140, year = 1993}
%
%  Notes:
%  ======
%
%  Note 1)
%
%  The original Jade algorithm/code deals with complex signals in
%  Gaussian noise white and exploits an underlying assumption that the
%  model of independent components actually holds.  This is a
%  reasonable assumption when dealing with some narrowband signals.
%  In this context, one may i) seriously consider dealing precisely
%  with the noise in the whitening process and ii) expect to use the
%  small number of significant eigenmatrices to efficiently summarize
%  all the 4th-order information.  All this is done in the JADE
%  algorithm.
%
%  In this implementation, we deal with real-valued signals and we do
%  NOT expect the ICA model to hold exactly.  Therefore, it is
%  pointless to try to deal precisely with the additive noise and it
%  is very unlikely that the cumulant tensor can be accurately
%  summarized by its first n eigen-matrices. Therefore, we consider
%  the joint diagonalization of the whole set of eigen-matrices.
%  However, in such a case, it is not necessary to compute the
%  eigenmatrices at all because one may equivalently use `parallel
%  slices' of the cumulant tensor.  This part (computing the
%  eigen-matrices) of the computation can be saved: it suffices to
%  jointly diagonalize a set of cumulant matrices.  Also, since we are
%  dealing with reals signals, it becomes easier to exploit the
%  symmetries of the cumulants to further reduce the number of
%  matrices to be diagonalized.  These considerations, together with
%  other cheap tricks lead to this version of JADE which is optimized
%  (again) to deal with real mixtures and to work `outside the model'.
%  As the original JADE algorithm, it works by minimizing a `good set'
%  of cumulants.
%
%  Note 2) 
%
%  The rows of the separating matrix B are resorted in such a way that
%  the columns of the corresponding mixing matrix A=pinv(B) are in
%  decreasing order of (Euclidian) norm.  This is a simple, `almost
%  canonical' way of fixing the indetermination of permutation.  It
%  has the effect that the first rows of the recovered signals (ie the
%  first rows of B*X) correspond to the most energetic *components*.
%  Recall however that the source signals in S=B*X have unit variance.
%  Therefore, when we say that the observations are unmixed in order
%  of decreasing energy, the energetic signature is found directly as
%  the norm of the columns of A=pinv(B).
%
%  Note 3) 
%  
%  In experiments where JADE is run as B=jadeR(X,m) with m varying in
%  range of values, it is nice to be able to test the stability of the
%  decomposition.  In order to help in such a test, the rows of B can
%  be sorted as described above. We have also decided to fix the sign
%  of each row in some arbitrary but fixed way.  The convention is
%  that the first element of each row of B is positive.
%
%
%  Note 4) 
%
%  Contrary to many other ICA algorithms, JADE (or least this version)
%  does not operate on the data themselves but on a statistic (the
%  full set of 4th order cumulant).  This is represented by the matrix
%  CM below, whose size grows as m^2 x m^2 where m is the number of
%  sources to be extracted (m could be much smaller than n).  As a
%  consequence, (this version of) JADE will probably choke on a
%  `large' number of sources.  Here `large' depends mainly on the
%  available memory and could be something like 40 or so.  One of
%  these days, I will prepare a version of JADE taking the `data'
%  option rather than the `statistic' option.
%
%


% JadeR.m ends here.

