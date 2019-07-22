clear;
clc;
close all;
data_dir = '/data/inprogress/kolagund/stromule_track/RawSplitChannels/segmented_old/';
lst = dir([data_dir,'chloroplast_*.tif']);
cmap1 = hsv(64);

parfor fno=1:numel(lst)
	fname = lst(fno).name;
	info = imfinfo([data_dir,fname]);
	save_filename = [data_dir,fname(1:end-4)];
	mkdir(save_filename);
	motion = {};
	disp(1);
	I1 = imread([data_dir,fname],1);
	s = regionprops(I1,'Centroid');
	pts = round(cat(1,s(:).Centroid));
	disp(size(pts))
	valid_inds=[];
	valid_inds{1} = true(size(pts,1),1);
	tracked_points=[];
	tracked_points{1} = pts;
	for i=1:numel(info)
		disp(i);
		I2 = imread([data_dir,fname],i);
		[D,~] = imregdemons(I1,I2,'DisplayWaitbar',false);
		motion{end+1} = D;
		I1 = I2;
		pts = max(1,round(pts + [D(sub2ind(size(D),pts(:,2),pts(:,1),ones(size(pts,1),1))), D(sub2ind(size(D),pts(:,2),pts(:,1),1+ones(size(pts,1),1)))]));
		pts(:,1) = min(pts(:,1),size(I2,2));
		pts(:,2) = min(pts(:,2),size(I2,1));
		
		disp(size(pts));
		valid_inds{end+1} = valid_inds{end} & I1(sub2ind(size(I1),pts(:,2),pts(:,1)));
		vi = valid_inds{end};
		It = I1;
		for j=1:size(pts,1)
			if vi(j)
				temp = regionprops(bwselect(I1,pts(j,1),pts(j,2)),'Centroid');
				pts(j,:) = round(temp(1).Centroid);
			end
		end
		tracked_points{end+1} = pts;
		
		% imshow(I1);
		% hold on;
		for j=1:size(pts,1)
			if vi(j)
				% plot(pts(j,1),pts(j,2),'*','color',uint8(255*cmap1(mod(j,64)+1,:)));
				dlmwrite([save_filename,'/',num2str(j),'.txt'],pts(j,:),'-append');
			end
		end
		% pause(0.1);
	end
end