%load 20news_w100

assigned = [1 2 3 4 1 4 5 8 5 3 6 5 6 4 1 6 6 1 6 ...
    4 6 7 6 4 8 8 2 6 5 6 6 2 3 8 6 8 1 1 2 2 4 5 6 ... % up to image
    8 3 3 3 7 8 2 7 6 7 1 6 7 7 6 7 2 5 4 7 1 6 5 2 5 8 5 ... % up to problem
    6 2 5 3 5 8 7 7 6 2 6 7 6 7 7 5 1 6 2 5 5 6 5 1 8 1 2 6 2 8];
z = assigned;

% for i = 1:100
%     if z(i) == 1
%         wordType = 'medical';
%     elseif z(i) == 2
%         wordType = 'sports';
%     elseif z(i) == 3
%         wordType = 'religion';
%     elseif z(i) == 4
%         wordType = 'cars';
%     elseif z(i) == 5
%         wordType = 'misc';
%     elseif z(i) == 6
%         wordType = 'computer';
%     elseif z(i) == 7
%         wordType = 'space';
%     elseif z(i) == 8
%         wordType = 'law/politics';
%     else
%         wordType = '0';
%     end
%     
%    fprintf('%s : %s\n',wordlist{i},wordType); 
% end