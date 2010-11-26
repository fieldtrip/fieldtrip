function []=tprod_testcases(testCases,debugin)
% This file contains lots of test-cases to test the performance of the tprod
% files vs. the matlab built-ins.
%
% 
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied
if ( nargin < 1 || isempty(testCases) ) testCases={'acc','timing','blksz'}; end;
if ( ~iscell(testCases) ) testCases={testCases}; end;
if ( nargin < 2 ) debugin=1; end;
global DEBUG; DEBUG=debugin;

% fprintf('------------------ memory execption test ------------\n');
% ans=tprod(randn(100000,1),1,randn(100000,1)',1); 
% ans=tprod(randn(100000,1),1,randn(100000,1)',1,'m');


%-----------------------------------------------------------------------
% Real , Real
if ( ~isempty(strmatch('acc',testCases)) ) 
fprintf('-------------------  Accuracy tests -------------------\n');

X=complex(randn(101,101,101),randn(101,101,101));
Y=complex(randn(101,101),randn(101,101));
fprintf('\n****************\n Double Real X, Double Real Y\n******************\n')
accuracyTests(real(X),real(Y),'dRdR');
fprintf('\n****************\n Double Real X, Double Complex Y\n******************\n')
accuracyTests(real(X),Y,'dRdC');
fprintf('\n****************\n Double Complex X, Double Real Y\n******************\n')
accuracyTests(X,real(Y),'dCdR');
fprintf('\n****************\n Double Complex X, Double Complex Y\n******************\n')
accuracyTests(X,Y,'dCdC');

fprintf('\n****************\n Double Real X, Single Real Y\n******************\n')
accuracyTests(real(X),single(real(Y)),'dRsR');
fprintf('\n****************\n Double Real X, Single Complex Y\n******************\n')
accuracyTests(real(X),single(Y),'dRsC');
fprintf('\n****************\n Double Complex X, Single Real Y\n******************\n')
accuracyTests(X,single(real(Y)),'dCsR');
fprintf('\n****************\n Double Complex X, Single Complex Y\n******************\n')
accuracyTests(X,single(Y),'dCsC');

fprintf('\n****************\n Single Real X, Double Real Y\n******************\n')
accuracyTests(single(real(X)),real(Y),'sRdR');
fprintf('\n****************\n Single Real X, Double Complex Y\n******************\n')
accuracyTests(single(real(X)),Y,'sRdC');
fprintf('\n****************\n Single Complex X, Double Real Y\n******************\n')
accuracyTests(single(X),real(Y),'sCdR');
fprintf('\n****************\n Single Complex X, Double Complex Y\n******************\n')
accuracyTests(single(X),Y,'sCdC');

fprintf('\n****************\n Single Real X, Single Real Y\n******************\n')
accuracyTests(single(real(X)),single(real(Y)),'sRsR');
fprintf('\n****************\n Single Real X, Single Complex Y\n******************\n')
accuracyTests(single(real(X)),single(Y),'sRsC');
fprintf('\n****************\n Single Complex X, Single Real Y\n******************\n')
accuracyTests(single(X),single(real(Y)),'sCsR');
fprintf('\n****************\n Single Complex X, Single Complex Y\n******************\n')
accuracyTests(single(X),single(Y),'sCsC');


fprintf('All tests passed\n');
end

% fprintf('-------------------  Timing tests -------------------\n');
if ( ~isempty(strmatch('timing',testCases)) ) 
%Timing tests
X=complex(randn(101,101,101),randn(101,101,101));
Y=complex(randn(101,101),randn(101,101));
fprintf('\n****************\n Real X, Real Y\n******************\n')
timingTests(real(X),real(Y),'RR');
fprintf('\n****************\n Real X, Complex Y\n******************\n')
timingTests(real(X),Y,'RC');
fprintf('\n****************\n Complex X, Real Y\n******************\n')
timingTests(X,real(Y),'CR');
fprintf('\n****************\n Complex X, Complex Y\n******************\n')
timingTests(X,Y,'CC');
end

% fprintf('------------------- Scaling  tests -------------------\n');
% scalingTests([32,64,128,256]);

if( ~isempty(strmatch('blksz',testCases)) )
fprintf('-------------------  Blksz tests -------------------\n');
blkSzTests([128 96 64 48 40 32 24 16 0],[128,256,512,1024,2048]);
end

return;


function []=accuracyTests(X,Y,str)
unitTest([str ' OuterProduct, [1],[2]'],tprod(X(:,1),1,Y(:,1),2,'m'),X(:,1)*Y(:,1).');
unitTest([str ' Inner product, [-1],[-1]'],tprod(X(:,1),-1,Y(:,1),-1,'m'),X(:,1).'*Y(:,1));
unitTest([str ' Matrix product, [1 -1],[-1 2]'],tprod(X(:,:,1),[1 -1],Y,[-1 2],'m'),X(:,:,1)*Y);
unitTest([str ' transposed matrix product, [-1 1],[-1 2]'],tprod(X(:,:,1),[-1 1],Y,[-1 2],'m'),X(:,:,1).'*Y);
unitTest([str ' Matrix frobenius norm, [-1 -2],[-1 -2]'],tprod(X(:,:,1),[-1 -2],Y,[-1 -2],'m'),sum(sum(X(:,:,1).*Y)));

unitTest([str ' transposed matrix frobenius norm, [-1 -2],[-2 -1]'],tprod(X(:,:,1),[-1 -2],Y,[-2 -1],'m'),sum(sum(X(:,:,1).'.*Y)));

unitTest([str ' ignored dims, [0 -2],[-2 2 1]'],tprod(Y(1,:),[0 -2],X(:,:,:),[-2 2 1],'m'),reshape(Y(1,:)*reshape(X,size(X,1),[]),size(X,2),size(X,3)).');


% Higher order matrix operations
unitTest([str ' spatio-temporal filter [-1 -2 1],[-1 -2]'],tprod(X,[-1 -2 1],Y,[-1 -2],'m'),reshape(Y(:).'*reshape(X,size(X,1)*size(X,2),size(X,3)),[size(X,3) 1]));

unitTest([str ' spatio-temporal filter (fallback) [-1 -2 1],[-1 -2]'],tprod(X,[-1 -2 1],Y,[-1 -2]),reshape(Y(:).'*reshape(X,size(X,1)*size(X,2),size(X,3)),[size(X,3) 1]));

unitTest([str ' spatio-temporal filter (order) [-1 -2],[-1 -2 1]'],tprod(Y,[-1 -2],X,[-1 -2 1]),reshape(Y(:).'*reshape(X,size(X,1)*size(X,2),size(X,3)),[size(X,3) 1]));

unitTest([str ' spatio-temporal filter (order) [-1 -2],[-1 -2 3]'],tprod(Y,[-1 -2],X,[-1 -2 3]),reshape(Y(:).'*reshape(X,size(X,1)*size(X,2),size(X,3)),[1 1 size(X,3)]));

unitTest([str ' transposed spatio-temporal filter [1 -2 -3],[-2 -3]'],tprod(X,[1 -2 -3],Y,[-2 -3],'m'),reshape(reshape(X,size(X,1),size(X,2)*size(X,3))*Y(:),[size(X,1) 1]));

unitTest([str ' matrix-vector product [-1 1 2][-1]'],tprod(X,[-1 1 2],Y(:,1),[-1],'m'),reshape(Y(:,1).'*reshape(X,[size(X,1) size(X,2)*size(X,3)]),[size(X,2) size(X,3)]));

unitTest([str ' spatial filter (fallback): [-1 2 3],[-1 1]'],tprod(X,[-1 2 3],Y,[-1 1]),reshape(Y.'*reshape(X,[size(X,1) size(X,2)*size(X,3)]),[size(Y,2) size(X,2) size(X,3)]));

unitTest([str ' spatial filter: [-1 2 3],[-1 1]'],tprod(X,[-1 2 3],Y,[-1 1],'m'),reshape(Y.'*reshape(X,[size(X,1) size(X,2)*size(X,3)]),[size(Y,2) size(X,2) size(X,3)]));

unitTest([str ' temporal filter [1 -2 3],[2 -2]'],tprod(X,[1 -2 3],Y(1,:),[2 -2],'m'),sum(X.*repmat(Y(1,:),[size(X,1) 1 size(X,3)]),2));

unitTest([str ' temporal filter [2 -2],[1 -2 3]'],tprod(Y(1,:),[2 -2],X,[1 -2 3],'m'),sum(X.*repmat(Y(1,:),[size(X,1) 1 size(X,3)]),2));

unitTest([str ' temporal filter [-2 2],[1 -2 3]'],tprod(Y(1,:).',[-2 2],X,[1 -2 3],'m'),sum(X.*repmat(Y(1,:),[size(X,1) 1 size(X,3)]),2));

unitTest([str ' temporal filter [-2 2],[1 -2 3]'],tprod(Y(1,:).',[-2 2],X,[1 -2 3],'m'),sum(X.*repmat(Y(1,:),[size(X,1) 1 size(X,3)]),2));

Xp=permute(X,[1 3 2]);
unitTest([str ' blk-code [-1 1 -2][-1 2 -2]'],tprod(X,[-1 1 -2],X,[-1 2 -2],'m'),reshape(Xp,[],size(Xp,3)).'*reshape(Xp,[],size(Xp,3)));

return;

function []=timingTests(X,Y,str);

% outer product simulation
fprintf([str ' OuterProduct [1][2]\n']);
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic; for i=1:1000; Z=tprod(X(:,1),1,Y(:,1),2); end
fprintf('%30s %gs\n','tprod',toc/1000); % = .05  / .01
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic; for i=1:1000; Zm=tprod(X(:,1),1,Y(:,1),2,'m'); end
fprintf('%30s %gs\n','tprod m',toc/1000); % = .05  / .01
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic, for i=1:1000;T=X(:,1)*Y(:,1).';end;
fprintf('%30s %gs\n','MATLAB',toc/1000); % = .03  / .01

% matrix product 
fprintf([str ' MatrixProduct [1 -1][-1 2]\n']);                    
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic, for i=1:1000;Z=tprod(X(:,:,1),[1 -1],Y,[-1 2]);end
fprintf('%30s %gs\n','tprod',toc/1000);% = .28 / .06
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic, for i=1:1000;Zm=tprod(X(:,:,1),[1 -1],Y,[-1 2],'m');end
fprintf('%30s %gs\n','tprod m',toc/1000);% = .28 / .06
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic, for i=1:1000;T=X(:,:,1)*Y;end
fprintf('%30s %gs\n','MATLAB',toc/1000);% = .17 / .06

% transposed matrix product simulation
fprintf([str ' transposed Matrix Product [-1 1][2 -1]\n']);
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic, for i=1:1000;Z=tprod(X(:,:,1),[-1 1],Y,[2 -1]);end
fprintf('%30s %gs\n','tprod',toc/1000);% =.3 / .06 
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic, for i=1:1000;Zm=tprod(X(:,:,1),[-1 1],Y,[2 -1],'m');end
fprintf('%30s %gs\n','tprod m',toc/1000);% =.3 / .06 
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic, for i=1:1000;T=X(:,:,1).'*Y.';end
fprintf('%30s %gs\n','MATLAB',toc/1000); % =.17  / .06

% Higher order matrix operations   % times: P3-m 1.8Ghz 2048k / P4 2.4Ghz 512k
fprintf([str ' spatio-temporal filter [-1 -2 1] [-1 -2]\n']);        
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:500;Z=tprod(X,[-1 -2 1],Y,[-1 -2]);end,
fprintf('%30s %gs\n','tprod',toc/500);% =.26 / .18
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:500;Zm=tprod(X,[-1 -2 1],Y,[-1 -2],'m');end,
fprintf('%30s %gs\n','tprod m',toc/500);% =.26 / .18
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:500;T=reshape(Y(:).'*reshape(X,size(X,1)*size(X,2),size(X,3)),[size(X,3) 1]);end,
fprintf('%30s %gs\n','MATLAB',toc/500); %=.21 / .18

fprintf([str ' transposed spatio-temporal filter [1 -2 -3] [-2 -3]\n']);    
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:50;Z=tprod(X,[1 -2 -3],Y,[-2 -3]);end,
fprintf('%30s %gs\n','tprod',toc/50);% =.27 / .28
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:50;Zm=tprod(X,[1 -2 -3],Y,[-2 -3],'m');end,
fprintf('%30s %gs\n','tprod m',toc/50);% =.27 / .28
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:50;T=reshape(reshape(X,size(X,1),size(X,2)*size(X,3))*Y(:),[size(X,1) 1]);end,
fprintf('%30s %gs\n','MATLAB',toc/50); %=.24 / .26

% MATRIX vector product
fprintf([str ' matrix-vector product [-1 1 2] [-1]\n']);
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:500; Z=tprod(X,[-1 1 2],Y(:,1),[-1]);end,
fprintf('%30s %gs\n','tprod',toc/500); %=.27 / .26
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:500; Zm=tprod(X,[-1 1 2],Y(:,1),[-1],'m');end,
fprintf('%30s %gs\n','tprod m',toc/500); %=.27 / .26
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:500; T=reshape(Y(:,1).'*reshape(X,[size(X,1) size(X,2)*size(X,3)]),[1 size(X,2) size(X,3)]);end,
fprintf('%30s %gs\n','MATLAB (reshape)',toc/500);%=.21 / .28
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:500; T=zeros([size(X,1),1,size(X,3)]);for k=1:size(X,3); T(:,:,k)=Y(1,:)*X(:,:,k); end,end,
fprintf('%30s %gs\n','MATLAB (loop)',toc/500); %=.49 /

% spatial filter  
fprintf([str ' Spatial filter: [-1 2 3],[-1 1]\n']);
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic;for i=1:500;Z=tprod(X,[-1 2 3],Y(:,1:2),[-1 1]);end;
fprintf('%30s %gs\n','tprod',toc/500);%=.39/.37
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic;for i=1:500;Zm=tprod(X,[-1 2 3],Y(:,1:2),[-1 1],'m');end;
fprintf('%30s %gs\n','tprod m',toc/500);%=.39/.37
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic;for i=1:500; T=reshape(Y(:,1:2).'*reshape(X,[size(X,1) size(X,2)*size(X,3)]),[2 size(X,2) size(X,3)]);end;
fprintf('%30s %gs\n','MATLAB',toc/500); 
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic;for i=1:500;T=zeros(2,size(X,2),size(X,3));for k=1:size(X,3);T(:,:,k)=Y(:,1:2).'*X(:,:,k);end;end;
fprintf('%30s %gs\n','MATLAB (loop)',toc/500); %=.76/.57

% temporal filter
fprintf([str ' Temporal filter: [1 -2 3],[2 -2]\n']);
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:50; Z=tprod(X,[1 -2 3],Y,[2 -2]);end
fprintf('%30s %gs\n','tprod',toc/50); %=.27 / .31
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:50; T=zeros([size(X,1),size(Y,1),size(X,3)]);for k=1:size(X,3); T(:,:,k)=X(:,:,k)*Y(:,:).'; end,end,
fprintf('%30s %gs\n','MATLAB (loop)',toc/50); %= .50 /
%tic,for i=1:50; T=zeros([size(X,1),size(Y,1),size(X,3)]); for k=1:size(Y,1); T(:,k,:)=sum(X.*repmat(Y(k,:),[size(X,1) 1 size(X,3)]),2);end,end;
%fprintf('%30s %gs\n','MATLAB (repmat)',toc/50); %=3.9 / 3.3

fprintf([str ' Temporal filter2: [1 -2 3],[-2 2]\n']);
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:50; Z=tprod(X,[1 -2 3],Y,[-2 2]);end
fprintf('%30s %gs\n','tprod',toc/50); %=.27 / .31
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:50; T=zeros([size(X,1),size(Y,1),size(X,3)]);for k=1:size(X,3); T(:,:,k)=X(:,:,k)*Y(:,:); end,end,
fprintf('%30s %gs\n','MATLAB (loop)',toc/50); %= .50 /

% Data Covariances                                           
fprintf([str ' Channel-covariance/trial(3) [1 -1 3] [2 -1 3]\n']); 
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic;for i=1:50;Z=tprod(X,[1 -1 3],[],[2 -1 3]);end;
fprintf('%30s %gs\n','tprod',toc/50);  %=8.36/7.6
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:50;T=zeros(size(X,1),size(X,1),size(X,3));for k=1:size(X,3); T(:,:,k)=X(:,:,k)*X(:,:,k)';end;end,
fprintf('%30s %gs\n','MATLAB (loop)',toc/50); %=9.66/7.0

% N.B. --- aligned over dim 1 takes 2x longer!                   
fprintf([str ' Channel-covariance/trial(1) [1 -1 2] [1 -1 3]\n']);
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic;for i=1:10;Z=tprod(X,[1 -1 2],[],[1 -1 3]);end;
fprintf('%30s %gs\n','tprod',toc/10); %=    /37.2
A=complex(randn(size(X)),randn(size(X))); % flush cache?
tic,for i=1:10;T=zeros(size(X,1),size(X,1),size(X,3));for k=1:size(X,3);T(k,:,:)=shiftdim(X(k,:,:))'*shiftdim(X(k,:,:));end;end,
fprintf('%30s %gs\n','MATLAB',toc/10);  %=17.2/25.8

return;


function []=scalingTests(Ns);
% Scaling test
fprintf('Simple test of the effect of the acc size\n');
for N=Ns;
fprintf('X=[%d x %d]\n',N,N*N);X=randn(N,N*N);Y=randn(N*N,1);
fprintf('[1 -1][-1]\n');
tic, for i=1:1000;Z=tprod(X,[1 -1],Y,[-1],'n');end
fprintf('%30s %gs\n','tprod',toc/1000);% = .28 / .06
tic, for i=1:1000;Z=tprod(X,[1 -1],Y,[-1],'mn');end
fprintf('%30s %gs\n','tprod m',toc/1000);% = .28 / .06
end

fprintf('Simple mat vs non mat tests\n');
for N=Ns;
fprintf('N=%d\n',N);X=randn(N,N,N);Y=randn(N,N);
fprintf('[1 -1 -2][-1 -2]\n');
tic, for i=1:1000;Z=tprod(X,[1 -1 -2],Y,[-1 -2]);end
fprintf('%30s %gs\n','tprod',toc/1000);% = .28 / .06
tic, for i=1:1000;Z=tprod(X,[1 -1 -2],Y,[-1 -2],'m');end
fprintf('%30s %gs\n','tprod m',toc/1000);% = .28 / .06
fprintf('[-1 -2 1][-1 -2]\n');
tic, for i=1:1000;Z=tprod(X,[-1 -2 1],Y,[-1 -2]);end
fprintf('%30s %gs\n','tprod',toc/1000);% = .28 / .06
tic, for i=1:1000;Z=tprod(X,[-1 -2 1],Y,[-1 -2],'m');end
fprintf('%30s %gs\n','tprod m',toc/1000);% = .28 / .06
fprintf('[-1 2 3][-1 1]\n');
tic, for i=1:100;Z=tprod(X,[-1 2 3],Y,[-1 1]);end
fprintf('%30s %gs\n','tprod',toc/100);% = .28 / .06
tic, for i=1:100;Z=tprod(X,[-1 2 3],Y,[-1 1],'m');end
fprintf('%30s %gs\n','tprod m',toc/100);% = .28 / .06
end

function []=blkSzTests(blkSzs,Ns);

fprintf('Blksz optimisation tests\n');
%blkSzs=[64 48 40 32 24 16 0];
tptime=zeros(numel(blkSzs+2));
for N=Ns;%[128,256,512,1024,2048];
  X=randn(N,N);Y=randn(size(X));
  for i=1:5;

     % use tprod without matlab
     for b=1:numel(blkSzs); blkSz=blkSzs(b);
      clear T;T=randn(1024,1024); % clear cache
      tic,for j=1:3;tprod(X,[1 -1],Y,[-1 2],'mn',blkSz);end;
		tptime(b)=tptime(b)+toc;
    end;

    % tprod with defaults & matlab if possible
    clear T;T=randn(1024,1024); % clear cache
    tic,for j=1:3;tprod(X,[1 -1],Y,[-1 2]);end;
	 tptime(numel(blkSzs)+1)=tptime(numel(blkSzs)+1)+toc;

    % the pure matlab code
    clear T;T=randn(1024,1024); % clear cache
    tic,for j=1:3; Z=X*Y;end;
	 tptime(numel(blkSzs)+2)=tptime(numel(blkSzs)+2)+toc;
  end;
  for j=1:numel(blkSzs);
    fprintf('blk=%d, N = %d -> %gs \n',blkSzs(j),N,tptime(j)); 
  end
  fprintf('tprod, N = %d -> %gs \n',N,tptime(numel(blkSzs)+1)); 
  fprintf('MATLAB, N = %d -> %gs \n',N,tptime(numel(blkSzs)+2));
  fprintf('\n');
end;
return;


% simple function to check the accuracy of a test and report the result
function [testRes,trueRes,diff]=unitTest(testStr,testRes,trueRes,tol)
global DEBUG;
if ( nargin < 4 ) 
   if ( isa(trueRes,'double') ) tol=1e-11; 
   elseif ( isa(trueRes,'single') ) tol=1e-5; 
   elseif ( isa(trueRes,'integer') ) 
      warning('Integer inputs!'); tol=1;       
   elseif ( isa(trueRes,'logical') ) tol=0;
   end
end
diff=abs(testRes-trueRes)./max(1,abs(testRes+trueRes));
fprintf('%60s = %0.3g ',testStr,max(diff(:)));
if ( max(diff(:)) > tol ) 
   if ( exist('mimage') )
      mimage(squeeze(testRes),squeeze(trueRes),squeeze(diff))
   end
   fprintf(' **FAILED***\n');   
   if ( DEBUG>0 ) 
      warning([testStr ': failed!']), 
      fprintf('Type return to continue\b');keyboard; 
   end;
else
   fprintf('Passed \n');
end
