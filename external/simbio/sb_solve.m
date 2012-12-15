function x = sb_solve(sysmat,vecb);

% SB_SOLVE
%
% $Id$

%scalen
disp('Scaling stiffnes matrix...')
dkond = 1./(sqrt(diag(sysmat)));
vecb = vecb.*dkond;
[indexi indexj s] = find(sysmat);
sys_size = size(sysmat,1);
clear sysmat;
s = (s.*dkond(indexi)).*dkond(indexj);
s(1) = 1;
disp('Preconditioning...')
L = sparse(indexi,indexj,s,sys_size,sys_size,length(s));
%partch
try
    L = ichol(L);
catch
    disp('Could not compute incomplete Cholesky-decompositon. Rescaling stiffness matrix...')
    alpha = 0.5d-6;
    alpha = alpha*8.d0;
    alpha = 1 / (alpha + 1);
    s = alpha*s;
    dia = find(indexi == indexj);
    s(dia) = (1./alpha)*s(dia);
    s(dia) = sqrt(s(dia));
    s(1) = 1;
    L = sparse(indexi,indexj,s,sys_size,sys_size,length(s));
    clear dia;
    L = ichol(L);
end
%startvektor
disp('Finding startvector...')
vecb_ = L \ (-vecb);
vecx = L' \ vecb_;
clear vecb_;
%sonstiges
sysmat = sparse(indexi,indexj,s,sys_size,sys_size,length(s));
sysmat = sysmat + sysmat' - sparse(1:sys_size,1:sys_size,diag(sysmat),sys_size,sys_size,sys_size);
clear indexi indexj s;
%dprod = sysmat * vecx;
%lösen
disp('Solving equation system...')
[x,fl,rr,it,rv]= pcg(sysmat,vecb,10e-9,5000,L,L',vecx);
%fl
rr
%it
%rescal
x = x.*dkond;
end
