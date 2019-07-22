function all_stromules = get_snake_pts(I)
cc = bwconncomp(I); %connected components

num_segments = 10; % number of parts/control points on the stromule
segments = 0:1/(num_segments-1):1;

all_stromules = cell(1,cc.NumObjects);
ccpixellist = cc.PixelIdxList; %list of pixels ids in each component
parfor i=1:cc.NumObjects
    [r,c] = ind2sub(size(I),ccpixellist{i}); %get (y,x) location from pixel ids
    [er,ec] = find(bwmorph(bwselect(I,c(1),r(1),8),'endpoints')); %find end points of the component
		%[br,bc] = find(bwmorph(bwselect(I,c(1),r(1),8),'branchpoints'));
		
	%if there are more than 2 end points find the ones that are most distant from each other
    D = pdist2([ec(:),er(:)],[ec(:),er(:)]);
    [p,q] = find(D==max(D(:)));
    
	%using end points initialize snake segments/parts/control points
    start_pt = [ec(p(1)), er(p(1))];
    end_pt = [ec(q(1)), er(q(1))];
    ln = end_pt - start_pt;
    ln_l = norm(ln);
    ln_d = ln/(eps+ln_l);
    
    x = start_pt(1) + ln_d(1)*segments*ln_l;
    y = start_pt(2) + ln_d(2)*segments*ln_l;
    
    if numel(r)>4
        all_stromules{i} = [x(:),y(:)];
    end
end

all_stromules = all_stromules(~cellfun('isempty',all_stromules)); 