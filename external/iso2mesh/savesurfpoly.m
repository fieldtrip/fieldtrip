function savesurfpoly(v,f,holelist,regionlist,p0,p1,fname,forcebox)
%
% savesurfpoly(v,f,holelist,regionlist,p0,p1,fname)
%
% save a set of surfaces into poly format (for tetgen)
%
% author: Qianqian Fang (fangq<at> nmr.mgh.harvard.edu)
% date: 2007/11/21
%
% input:
%      v: input, surface node list, dimension (nn,3)
%         if v has 4 columns, the last column specifies mesh density near each node
%      f: input, surface face element list, dimension (be,3)
%      holelist: list of holes, each hole is represented by an internal point
%      regionlist: list of regions, similar to holelist
%      p0: coordinate of one of the end of the bounding box
%      p1: coordinate for the other end of the bounding box
%      fname: output file name
%      forcebox: non-empty: add bounding box, []: automatic
%                if forcebox is a 8x1 vector, it will be used to 
%                specify max-edge size near the bounding box corners
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
dobbx=0;
if(nargin>=8)
	dobbx=~isempty(forcebox) & all(forcebox);
end

if(~iscell(f) & size(f,2)==4)
    faceid=f(:,4);
    f=f(:,1:3);
end

if(~iscell(f))
    edges=surfedge(f);
else
    edges=[];
end
bbxnum=0;

nodesize=[];
if(size(v,2)==4)
   nodesize=v(:,4);
   v=v(:,1:3);
end
node=v;
loopid=[];
if(~isempty(edges))
    loops=extractloops(edges);
    if(length(loops)<3)
        error('degenerated loops detected');
    end
    seg=[0,find(isnan(loops))];
    segnum=length(seg)-1;
    newloops=[];
    for i=1:segnum
       if(seg(i+1)-(seg(i)+1)==0) continue; end
       newloops=[newloops nan bbxflatsegment(node,loops(seg(i)+1:seg(i+1)-1))];
    end
    loops=[newloops nan];

    seg=[0,find(isnan(loops))];
    segnum=length(seg)-1;
    bbxnum=6;
    loopcount=zeros(bbxnum,1);
    loopid=zeros(segnum,1);
    for i=1:segnum     % walk through the edge loops
        subloop=loops(seg(i)+1:seg(i+1)-1);
        if(isempty(subloop)) continue; end
        boxfacet=find(sum(abs(diff(v(subloop,:))))<1e-8); % find a flat loop
        if(length(boxfacet)==1)   % if the loop is flat along x/y/z dir
            bf=boxfacet(1);    % no degeneracy allowed
            if(sum(abs(v(subloop(1),bf)-p0(bf)))<1e-2)
                loopcount(bf)=loopcount(bf)+1;
                v(subloop,bf)=p0(bf);
                loopid(i)=bf;
            elseif(sum(abs(v(subloop(1),bf)-p1(bf)))<1e-2)
                loopcount(bf+3)=loopcount(bf+3)+1;
                v(subloop,bf)=p1(bf);
                loopid(i)=bf+3;
            end
        end
    end
end

if(dobbx & isempty(edges))
    bbxnum=6;
    loopcount=zeros(bbxnum,1);	
end

if(dobbx|~isempty(edges))
    nn=size(v,1);
    boxnode=[p0;p1(1),p0(2:3);p1(1:2),p0(3);p0(1),p1(2),p0(3);
              p0(1:2),p1(3);p1(1),p0(2),p1(3);p1;p0(1),p1(2:3)];
    boxelem=[
        4 nn nn+3 nn+7 nn+4;   % x=xmin
        4 nn nn+1 nn+5 nn+4;   % y=ymin
        4 nn nn+1 nn+2 nn+3;   % z=zmin
        4 nn+1 nn+2 nn+6 nn+5; % x=xmax
        4 nn+2 nn+3 nn+7 nn+6; % y=ymax
        4 nn+4 nn+5 nn+6 nn+7];% z=zmax

    node=[v;boxnode];
end

node=[(0:size(node,1)-1)',node];

fp=fopen(fname,'wt');
fprintf(fp,'#node list\n%d 3 0 0\n',length(node));
fprintf(fp,'%d %f %f %f\n',node');

if(~iscell(f))
    fprintf(fp,'#facet list\n%d 1\n',length(f)+bbxnum);
    elem=[3*ones(length(f),1),f-1,ones(length(f),1)];
    if(~isempty(elem))
	fprintf(fp,'1 0\n%d %d %d %d %d\n',elem');
    end
else % if the surface is recorded as a cell array
    totalplc=0;
    for i=1:length(f)
        if(~iscell(f{i}))
            totalplc=totalplc+size(f{i},1);
        else
            totalplc=totalplc+size(f{i}{1},1);
        end
    end
    fprintf(fp,'#facet list\n%d 1\n',totalplc+bbxnum);
    for i=1:length(f)
        plcs=f{i};
        faceid=-1;
        if(iscell(plcs)) % if each face is a cell, use plc{2} for face id
            if(length(plcs)>1)
                faceid=plcs{2};
            end
            plcs=plcs{1};
        end
        for row=1:size(plcs,1);
         plc=plcs(row,:);
         if(any(isnan(plc))) % we use nan to separate outter contours and holes
            holeid=find(isnan(plc));
            if(faceid>0) 
                fprintf(fp,'%d %d %d\n%d',length(holeid)+1,length(holeid),faceid,holeid(1)-1);
            else
                fprintf(fp,'%d %d\n%d',length(holeid)+1,length(holeid),holeid(1)-1);
            end
            fprintf(fp,'\t%d',plc(1:holeid(1)-1)-1);
            fprintf(fp,'\t1\n');
            for j=1:length(holeid)
                if(j==length(holeid))
                    fprintf(fp,'%d',length(plc(holeid(j)+1:end)));
                fprintf(fp,'\t%d',plc(holeid(j)+1:end)-1);
                else
                    fprintf(fp,'%d',length(plc(holeid(j)+1:holeid(j+1)-1)));
                    fprintf(fp,'\t%d',plc(holeid(j)+1:holeid(j+1)-1)-1);
                end
                fprintf(fp,'\t1\n');
            end
            for j=1:length(holeid)
                if(j==length(holeid))
                    fprintf(fp,'%d %f %f %f\n',j,mean(node(plc(holeid(j)+1:end),2:4)));
                else
                    fprintf(fp,'%d %f %f %f\n',j,mean(node(plc(holeid(j)+1:holeid(j+1)-1),2:4)));
                end
            end
         else
	    	if(faceid>0)
                fprintf(fp,'1 0 %d\n%d',faceid,length(plc));
            else
                fprintf(fp,'1 0\n%d',length(plc));
            end
    		fprintf(fp,'\t%d',plc-1);
            fprintf(fp,'\t1\n');
         end
        end
    end
end

if(dobbx|~isempty(edges))
    for i=1:bbxnum
        fprintf(fp,'%d %d 1\n',1+loopcount(i),loopcount(i));
        fprintf(fp,'%d %d %d %d %d\n',boxelem(i,:));
        if(~isempty(edges) & loopcount(i) &~isempty(find(loopid==i)))
            endid=find(loopid==i);
            for k=1:length(endid)
                j=endid(k);
                subloop=loops(seg(j)+1:seg(j+1)-1);
                fprintf(fp,'%d ',length(subloop));
                fprintf(fp,'%d ',subloop-1);
                fprintf(fp,'\n');
            end
            for k=1:length(endid)
                j=endid(k);
                subloop=loops(seg(j)+1:seg(j+1)-1);
                fprintf(fp,'%d %f %f %f\n',k,internalpoint(v,subloop)); %mean(v(subloop,:)));
            end
        end
    end
end

if(size(holelist,1))
        fprintf(fp,'#hole list\n%d\n',size(holelist,1));
        for i=1:size(holelist,1)
                fprintf(fp,'%d %f %f %f\n', i, holelist(i,:));
        end
else
	fprintf(fp,'#hole list\n0\n');
end

if(size(regionlist,1))
	fprintf(fp,'#region list\n%d\n',size(regionlist,1));
	for i=1:size(regionlist,1)
		fprintf(fp,'%d %f %f %f %d\n', i, regionlist(i,:),i);
	end
end
fclose(fp);

if(~isempty(nodesize))
	if(size(nodesize,1)+size(forcebox(:),1)==size(node,1))
		nodesize=[nodesize;forcebox(:)];
	end
	fid=fopen(regexprep(fname,'\.poly$','.mtr'),'wt');
	fprintf(fid,'%d 1\n',size(nodesize,1));
	fprintf(fid,'%f\n',nodesize);
	fclose(fid);
end
