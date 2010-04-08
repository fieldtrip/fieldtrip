function str = avgoverlabel(label)

str = sprintf('%s, ', label{:});
str = str(1:(end-2));
str = sprintf('mean(%s)', str);
str = {str};
