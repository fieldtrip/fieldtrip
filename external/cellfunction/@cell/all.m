function bool = all(x)

s = size2(x,[],'cell');
n = sum(cumprod(s,2));
x = x~=0;

bool = sum(sum(x))==n(end);
