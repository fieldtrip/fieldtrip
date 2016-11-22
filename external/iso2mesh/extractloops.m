function loops=extractloops(edges)
%
% loops=extractloops(edges)
%
% extract individual loop or polyline segment from a collection of edges
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2007/11/21
%
% input:   
%    edges:  two column matrix recording the starting/ending 
%             points of all edge segments
%
% output:
%    loops:  output, a single vector separated by NaN, each segment
%             is a 3D polyline or loop consisted of node IDs
%
% example:
%    edges=[1 2;2 3;1 4;3 4;7 3;1 9;5 6;6 7;10 9; 8 10;1 8;9 3;11 11;11 12];
%    loops=extractloops(edges)
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

loops=[];
loops=[loops,edges(1,:)];
loophead=edges(1,1);
loopend=edges(1,end);
edges(1,:)=[];

while(~isempty(edges))
    idx=[find(edges(:,1)==loopend)',find(edges(:,2)==loopend)'];
    if(length(idx)>1) % when a node with multiple connection found
        idx=idx(1);   % take the first connection and continue
    end
    if(isempty(idx)) % when an open-line segment gets to one end
        % when both open ends are found
        if(isempty([find(edges(:,1)==loophead)',find(edges(:,2)==loophead)']))
            loops=[loops,nan];
            loops=[loops,edges(1,:)];
            loophead=edges(1,1);
            loopend=edges(1,end);
            edges(1,:)=[];
        else % only the first open end is found, flip and trace the other
            [loophead, loopend]=deal(loopend, loophead);
            lp=fliplr(loops);
            seg=find(isnan(lp),1);
            if(isempty(seg))
                loops=lp;
            else
                loops=[loops(1:end-seg(1)+1) lp(1:seg(1)-1)];
            end
        end
        continue;    
    end
    if(length(idx)==1) % tracing along a single line thread
        idx=idx(1);
        newend=setdiff(edges(idx,:),loopend);
        if(newend==loophead)  % when a loop is found
            loops=[loops loophead nan];
            edges(idx,:)=[];
            if(size(edges,1)==0) break; end
            loops=[loops,edges(1,:)];
            loophead=edges(1,1);
            loopend=edges(1,end);
            edges(1,:)=[];
            continue;
        else
            loops=[loops,newend];
        end
        loopend=newend;
        edges(idx,:)=[];
    end
end
    
