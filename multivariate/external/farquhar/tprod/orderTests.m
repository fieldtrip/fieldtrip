function [Z2]=orderTests(X,Y,sz)
if ( nargin < 1 ) X=[]; end;
if( nargin<2 || isempty(Y) ) if(~isstr(X)) Y=randn(size(X)); else Y=[]; end;end
if ( nargin < 3 || isempty(sz) ) sz=[11 29 51]; end;
if ( isempty(X)||isstr(X) || isempty(Y)||isstr(Y) ) 
   if (isstr(X)) X=feval(X,randn(sz)); else X=randn(sz); end;
   if (isstr(Y)) Y=feval(Y,randn(sz)); else Y=randn(sz); end;
end
if ( isa(X,'single') ) cX='s'; else cX='d'; end;
if ( isa(Y,'single') ) cY='s'; else cY='d'; end;

% linear over last dim, ip over 1st, op over 2nd
Z=[];for i=1:size(X,3); Z(:,:,i)=X(:,:,i).'*Y(:,:,i); end;
Z2=tprod(X,[-1 1 3],Y,[-1 2 3],'m');
unitTest([cX 'X,[-1 1 3], ' cY 'Y,[-1 2 3]'],Z2,Z,1e-4);

Z=[];for i=1:size(X,3); Z(:,:,i)=Y(:,:,i).'*X(:,:,i); end;
Z2=tprod(X,[-1 2 3],Y,[-1 1 3],'m');
unitTest([cX 'X,[-1 2 3], ' cY 'Y,[-1 1 3]'],Z2,Z,1e-4);


% linear over last dim, op over 1st, ip over 2nd
Z=[];for i=1:size(X,3); Z(:,:,i)=(Y(:,:,i)*X(:,:,i).').'; end;
Z2=tprod(X,[1 -2 3],Y,[2 -2 3],'m');
unitTest([cX 'X,[1 -2 3], ' cY 'Y,[2 -2 3]'],Z2,Z,1e-4);

Z=[];for i=1:size(X,3); Z(:,:,i)=X(:,:,i)*Y(:,:,i).'; end;
z2=tprod(X,[2 -2 3],Y,[1 -2 3],'m');
unitTest([cX 'X,[2 -2 3], ' cY 'Y,[1 -2 3]'],Z2,Z,1e-4);

return;

% linear over last dim, op over 1+2nd, op over 2+1
Z2=tprod(permute(X,[2 1 3]),[-1 2 3],        Y,         [1 -1 3]);
unitTest('X,[-1 2 3],Y,[1 -1 3]',Z2,Z);

Z2=tprod(        X,         [1 -1 3],permute(Y,[2 1 3]),[-1 2 3]);
unitTest('X,[1 -1 3],Y,[-1 2 3]',Z2,Z)

%----------------------------------------

szs={[11 29 51] [29 11 13] [51 277 93]};
Xs=cell(20,numel(szs)); Ys=cell(20,numel(szs)); Zs=cell(20,numel(szs));
for i=1:1000;
   szi=ceil(rand(1,1)*numel(szs));
   tri=ceil(rand(1,1)*size(Xs,1));
   if ( ~isempty(Xs{tri,szi}) ) % use to check validity
      sum(Xs{tri,szi}(:));sum(Xs{tri,szi}(:));sum(Zs{tri,szi}(:));
      Zs{tri,szi}=orderTests(Xs{tri,szi},Ys{tri,szi});
   end
   % Overwrite with new
   Xs{tri,szi}=randn(szs{szi});Ys{tri,szi}=randn(szs{szi});
   if ( randn(1,1)>0 ) Xs{tri,szi}=single(Xs{tri,szi}); end
   if ( randn(1,1)>0 ) Ys{tri,szi}=single(Ys{tri,szi}); end
end
