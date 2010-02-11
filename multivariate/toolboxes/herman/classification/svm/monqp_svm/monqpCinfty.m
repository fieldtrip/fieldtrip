% This function has been adopted from S. Canu's work
% (scanu@insa-rouen.fr,INSA de Rennes) - under GNU license

function [xnew, lambda, ind] = monqpCinfty(H,c,A,b,l,verbose,X,ps,xinit)
% 
% min 1/2  x' H x - c' x                
%  x 
% contraints:   A' x = b   
%               y'*x = 0 
%               0 <= x 
%
 
 
%-------------------------------------------------------------------------- 
%                                verifications 
%-------------------------------------------------------------------------- 
[n,d] = size(H); 
[nl,nc] = size(c); 
[nlc,ncc] = size(A); 
[nlb,ncb] = size(b); 
if d ~= n 
	error('H must be a squre matrix n by n'); 
end 
if nl ~= n 
	error('H and c must have the same number of row'); 
end 
 
if nlc ~= n 
	error('H and A must have the same number of row'); 
end 
if nc ~= 1 
	error('c must be a row vector'); 
end 
if ncb ~= 1 
	error('b must be a row vector'); 
end 
if ncc ~= nlb  
	error('A'' and b  must have the same number of row'); 
end 
 
if nargin < 5		% default value for the regularisation parameter 
	l = 0;				 
end; 
 
if nargin < 6		% default value for the display parameter 
	verbose = 0;				 
end; 

if exist('xinit','var')~=1
    xinit=[];
end;
 
 
fid = 1; %default value, curent matlab window 
%-------------------------------------------------------------------------- 
 

%-------------------------------------------------------------------------- 
%                       I N I T I A L I S A T I O N  
%-------------------------------------------------------------------------- 
 
OO = zeros(ncc);
H = H+l*eye(length(H)); % preconditionning

xnew = -1;ness = 0;
if isempty(xinit)
while (min(xnew) < 0 ) & ness < 100
   ind = randperm(n);
   ind = sort(ind(1:ncc)');
   aux=[H(ind,ind) A(ind,:);A(ind,:)' OO];
   aux= aux+l*eye(size(aux));
   newsol = aux\[c(ind) ; b];
   xnew = newsol(1:length(ind));
   ness = ness+1;
end;
x = xnew;
lambda = newsol(length(ind)+1:length(ind)+ncc);
else  % start with a predefined x

%     ind=find(xinit);
%        
%    aux=[H(ind,ind) A(ind,:);A(ind,:)' OO];
%     aux= aux+l*eye(size(aux));
%     newsol = aux\[c(ind) ; b];
%     xnew = newsol(1:length(ind));
%     lambda = newsol(length(ind)+1:length(ind)+ncc);
%     
%     x = xnew;
    
 
    ind=find(xinit>0);    
    x = xinit(ind);
    lambda=ones(ncc,1);
    
    
end;


%ind = 1;
%x = C/2; xnew = x;
%keyboard;
Hini = H;Aini=A;cini=c;bini=b;

eyen = l*eye(n);

%-------------------------------------------------------------------------- 
%                            M A I N   L O O P 
%-------------------------------------------------------------------------- 

Jold = 10000000000000000000;  
C = Jold; % for the cost function
if verbose ~= 0 
  disp('      Cost     Delta Cost  #support  #up saturate'); 
  nbverbose = 0; 
end 

STOP = 0;
while STOP ~= 1



indd = zeros(n,1);indd(ind)=ind;nsup=length(ind);indsuptot=[];
[J,yx] = cout(H,x,b,C,indd,c,A,b,lambda); 
if verbose ~= 0 
  nbverbose = nbverbose+1; 
  if nbverbose == 20 
  disp('      Cost     Delta Cost  #support  #up saturate'); 
  nbverbose = 0; 
end 
if Jold == 0 
   fprintf(fid,'| %11.4e | %8.4f | %6.0f | %6.0f |\n',[J (Jold-J) nsup length(indsuptot)]); 
   elseif  (Jold-J)>0
fprintf(fid,'| %11.4e | %8.4f | %6.0f | %6.0f |\n',[J min((Jold-J)/abs(Jold),99.9999) nsup length(indsuptot)]); 
   else
fprintf(fid,'| %11.4e | %8.4f | %6.0f | %6.0f | bad mouve \n',[J max((Jold-J)/abs(Jold),-99.9999) nsup length(indsuptot)]); 
end 
end 
Jold = J; 
 
% nouvele resolution du système

% newsol = [H(ind,ind)+l*eye(length(ind)) A(ind,:);A(ind,:)' OO]\[c(ind) ; b];
% xnew = newsol(1:length(ind));
% lambda = newsol(length(ind)+1:length(newsol));

% %keyboard
% U = chol(H(ind,ind));      % CHOLESKI INSIDE!
% At=A(ind,:)';
% Ae=A(ind,:);
% ce = c(ind);
% Ut = U';
% M = At*(U\(Ut\Ae));
% d = U\(Ut\ce);
% d = (At*d - b);    % second membre  (homage au gars OlgaZZZZZZ qui nous a bien aidé)
%                    % On appelle le gc fabuleux de mister OlgaZZZ Hoduc
%                    % lambdastart = rand([m,1]);
%                    % [lambda] = gradconj(lambdastart,M*M,M*d,MaxIterZZZ,tol,1000);
% lambda = M\d;
%                    % On reconstitue le vecteur en résolvant le système linéaire Hx = c-Alambda
% xnew = U\(Ut\(ce-Ae*lambda));
    auxH=H(ind,ind);
    [U,testchol] = chol(auxH);
    At=A(ind,:)';
    Ae=A(ind,:);
    ce = c(ind);
    if testchol==0
        Ut = U';
        M = At*(U\(Ut\Ae));
        d = U\(Ut\ce);
        d = (At*d - b);   % second membre  (homage au gars OlgaZZZZZZ qui nous a bien aidé)
                            % On appelle le gc fabuleux de mister OlgaZZZ Hoduc
                            % lambdastart = rand([m,1]);
                            % [lambda] = gradconj(lambdastart,M*M,M*d,MaxIterZZZ,tol,1000);
        lambda = M\d;
        % On reconstitue le vecteur en résolvant le système linéaire Hx = c-Alambda
        xnew = U\(Ut\(ce-Ae*lambda));
    else
        
        % if non definite positive;
        
        
        auxM=At*(auxH\Ae);
        M=auxM'*auxM;
        d=auxH\ce;
        d=At*d-b;
        lambda=M\d;
        xnew=auxH\(ce-Ae*lambda);
    end;



%-------------------------------------------------------------------------------------------------------

[minxnew minpos]=min(xnew);

if min(xnew) < 0  &  length(ind) > ncc
                                  % on projette dans le domaine admissible et on sature la contrainte associée
%                                   -------------------------------------------------------------------------
     d = xnew - x;                % projection direction
     indad = find(xnew < 0);
     [t indmin] = min(-x(indad)./d(indad));

     x = x + t*d;                 % projection into the admissible set

     ind(indad(indmin)) = [];            % the associate variable is set to 0
     x(indad(indmin)) = [];              % the associate variable is set to 0

else
     xt = zeros(n,1);             % keyboard;
     xt(ind) = xnew;              % keyboard;

     mu = H*xt - c + A*lambda;    % calcul des multiplicateurs de lagrange associées aux contraintes

     indsat = 1:n;                           % on ne regarde que les contraintes saturées
     indsat(ind) = [];

     [mm mpos]=min(mu(indsat));

     if mm < -sqrt(eps)            % on enleve une contrainte active (saturée) 
         ind = sort([ind ; indsat(mpos)]); % et on rajoute en inactif
         x = xt(ind);
     else
         STOP = 1;                % AILLLET !
     end

end
end


%-------------------------------------------------------------------------- 
%                        E N D   M A I N   L O O P 
%-------------------------------------------------------------------------- 
 
