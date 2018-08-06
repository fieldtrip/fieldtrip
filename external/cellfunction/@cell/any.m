function bool = any(x)

s = size2(x,[],'cell');
n = sum(cumprod(s,2));
x = x==1;

bool = sum(sum(x))==n(end);
