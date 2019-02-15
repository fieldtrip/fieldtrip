function R=cpd_R(a,b,g)

if nargin==1
    R=rot(a);
end

if nargin==3

    R1=eye(3);
    R2=eye(3);
    R3=eye(3);

    R1(1:2,1:2)=rot(a);
    R2([1 3],[1 3])=rot(b);
    R3(2:3,2:3)=rot(g);

    R=R1*R2*R3;

end



function R=rot(a)

ca=cos(a);
sa=sin(a);

R=[ca -sa;
    sa ca];