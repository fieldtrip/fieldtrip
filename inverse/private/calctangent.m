function [tanu, tanv] = calctangent(RDip);

%% Based on calcrads.m, only difference is that RDip is alread
%% with respect to the sphere origin in calctangent.m
%% MODIFIED 13th JAN 2005 MATT BROOKES

x=RDip(1);
y=RDip(2);
z=RDip(3);
r=sqrt(x*x+y*y+z*z);

if (x==0) & (y==0)
  tanu(1)=1.0; tanu(2)=0; tanu(3)=0;
  tanv(1)=0; tanv(2)=1.0; tanv(3)=0;
else
  RZXY= -(r-z)*x*y;
  X2Y2= 1/(x*x+y*y);
  
  tanu(1)= (z*x*x + r*y*y) * X2Y2/r;
  tanu(2)= RZXY * X2Y2/r;
  tanu(3)= -x/r;
  
  tanv(1)= RZXY * X2Y2/r;
  tanv(2)= (z*y*y + r*x*x) * X2Y2/r;
  tanv(3)= -y/r;
end;  

  
