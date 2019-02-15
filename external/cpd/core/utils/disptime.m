% Function for convinient display of time from 'toc'
% use disptime(toc)
%
% Andriy Myronenko
% Feb 10, 2006

function disptime(t)

if t<60
    s=t;
    m=0;h=0;
end

if (t>=60) && (t<3600)
    m=floor(t/60);
    s=t-m*60;
    h=0;
end

if t>=3600
    h=floor(t/3600);
    m=floor((t-3600*h)/60);
    s=t-h*3600-m*60;
end

disp(['Time =' num2str(h) ' hours ' num2str(m) ' minutes ' num2str(s) ' seconds ']);