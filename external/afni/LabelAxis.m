function LabelAxis(ha, xlbls, ylbls)
%Puts text labels on axis of graph
%example:
%figure(1); plot ([1 2], [4 5]);
%LabelAxis(gca,{'Low' 'Hi'}, {'four' 'five'})
 
if (~isempty(xlbls) & ~iscellstr(xlbls)),
   fprintf(2,'Need array of cellstrings for x input');
end
if (~isempty(ylbls) & ~iscellstr(ylbls)),
   fprintf(2,'Need array of cellstrings for y input');
end

lx = length(xlbls);
ly = length(ylbls);

if (lx > 0),
   set(ha,'XTick',[1:1:lx]);
   set(ha,'XTickLabel',char(xlbls));
end
if (ly > 0),
   set(ha,'YTick',[1:1:ly]);
   set(ha,'YTickLabel',char(ylbls));
end

return;
