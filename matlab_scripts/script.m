clear;
close all;
filename = '1.gif';
for imno=0:59
    I = uint8(imread(['1/t',num2str(imno),'.png']));
    I_s = bwmorph(im2bw(I),'skel',Inf);
    RI = imread('1/1.lsm',imno*2+1);
    RI_s = repmat(RI(:,:,2),1,1,3);
    cc = bwconncomp(I_s);
    
    num_segments = 10;
    segments = 0:1/(num_segments-1):1;
    
    spl_pts = 0:0.001:1;

    for i=1:cc.NumObjects
        [r,c] = ind2sub(size(I),cc.PixelIdxList{i});
        [er,ec] = find(bwmorph(bwselect(I_s,c(1),r(1),8),'endpoints'));
        D = pdist2([ec(:),er(:)],[ec(:),er(:)]);
        [p,q] = find(D==max(D(:)));
        
        start_pt = [ec(p(1)), er(p(1))];
        end_pt = [ec(q(1)), er(q(1))];
        ln = end_pt - start_pt;
        ln_l = norm(ln);
        ln_d = ln/(eps+ln_l);
        
        [prj,idxs] = sort((( [c,r] - repmat(start_pt,numel(r),1) )*ln_d')/ln_l);
        idxs(prj>1)=[];

        if numel(r)>4
            inds = round(1:(numel(r)-1)/(num_segments-1):numel(r));
            x = spline(segments,c(inds),spl_pts);
            y = spline(segments,r(inds),spl_pts);
            RI_s = insertMarker(RI_s,[c(inds),r(inds)],'size',2);
            RI_s = insertShape(RI_s,'Line',reshape([c(inds)';r(inds)'],1,[]),'color','red');
        end
    end
    [RI_s,map] = rgb2ind(RI_s,256);
    if imno == 0 
          imwrite(RI_s,map,filename,'gif', 'Loopcount',inf); 
    else 
          imwrite(RI_s,map,filename,'gif','WriteMode','append'); 
    end 
end