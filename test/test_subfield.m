function test_subfield

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY issubfield getsubfield setsubfield

a.b.c = 1;
assert(issubfield(a,'b.c')==true);
assert(issubfield(a,'b')==true);
assert(issubfield(a.b,'c')==true);
assert(issubfield(a,'b.d')==false);
assert(issubfield(a.b,'d')==false);
assert(issubfield(a,'d.c')==false);
assert(issubfield(a,'d')==false);

assert(getsubfield(a.b,'c')==1);
assert(getsubfield(a,'b.c')==1);

a2 = setsubfield(a,'b.c',2);
assert(a2.b.c==2);
a2.b = setsubfield(a2.b,'c',4);
assert(a2.b.c==4);

