clear;
path = '/raid1/inprogress/kolagund/stromule_track/RawSplitChannels/';
lst = dir([path, '*.tif']);
for fno = 1:numel(lst)
    fname = lst(fno).name;
    disp([fname,'...']);
    info = imfinfo([path,fname]);
    [xx,yy] = meshgrid(1:info(1).Width,1:info(1).Height);
    for i=1:3:numel(info) %num time frames
				disp((i-1)/3)
        %%read image
        %if ~isempty(strfind(fname,'311'))						
				Ii = bfilter2(im2double(imread([path,fname],i)));
				Io = bfilter2(im2double(imread([path,fname],i+1)));
        %else
        %    Ii = im2double(imread([path,fname],i));
        %    Io = im2double(imread([path,fname],i+1));
        %end
        %Ti = Io;
        
        Ii = min(1,Ii/(0.9*max(Ii(:))));
				%Io = Io./Ii;
        Io = min(1,Io/(0.9*max(Io(:))));
        
        %%detect chloroplast: identify initial mask by thresholding red
        %%channel
        
        %%(thresh1), find regions connected to initial mask in cyan channel
        %%that are brighter than thresh2
        thresh1 = graythresh(Ii);
        thresh2 = graythresh(Io);
				bw1 = bwareaopen(im2bw(Ii,thresh1),15);
        [r,c] = find(bw1);
        mask = bwselect(im2bw(Io,thresh2),c,r,8);
        bw = imfill(imopen(mask,ones(11)),'holes');
        
        %%watershed to separate connected chloroplasts:
        %%http://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/?s_tid=srchtitle
				
				%D = -bwdist(~bw1);
        %imask = imextendedmin(-D,1);
        %D = imimposemin(D,imask);
        %Ld = watershed(D);
        %bw1(Ld == 0) = 0;
				%D = -bwdist(~bw1);
        %imask = imextendedmin(-D,1);
				
				%props = regionprops(bw1,'Centroid')
				%centroid = round(cat(2,props.Centroid));
				D = -bwdist(~bw);
        imask = imextendedmin(D,1);
        D = imimposemin(D,imask);
        Ld = watershed(D);
        bw(Ld == 0) = 0;
        Cmask = bwareafilt(bw,[200,2000]);
				
				if i==1
           imwrite(Cmask,['chloroplast_',fname]);
           %imwrite([mmask;nnmask],[path,'segment_conn_',fname,'.tif']);
        else
           imwrite(Cmask,['chloroplast_',fname],'writemode','append');
           %imwrite([mmask;nnmask],[path,'segment_conn_',fname,'.tif'],'writemode','append');
        end
		end
end	