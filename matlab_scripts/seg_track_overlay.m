blank =  uint8(zeros(600,1024));

%img33 = gifread('YY13_p50-HA_01_Maximumintensityprojection.gif');

data_dir = '/data/inprogress/kolagund/stromule_track/';

sample = 'YY13_p50-HA_04_Maximumintensityprojection';

tracks = dir([data_dir,sample ,'_new/updated/*.txt']);
    
strom_img = repmat(blank,1,1,3,60);
strom_tip = repmat(blank,1,1,3);
for tno=1:numel(tracks)
    track = tracks(tno).name;
    st = dlmread([data_dir, sample, '_new/updated/', track]);
    xs = st(:,1:10);%xs = st(1:3:end,:);
    ys = st(:,11:20);%ys = st(2:3:end,:);
    zs = st(:,21:end);
    for j = 1:size(zs,1)
        strom_img(:,:,:,zs(j,1)+1) = insertShape(strom_img(:,:,:,zs(j,1)+1),'Line',reshape([xs(j,:);ys(j,:)],1,[]),'LineWidth',3,'color',[0,255,0]);
        strom_img(:,:,:,zs(j,1)+1) = insertMarker(strom_img(:,:,:,zs(j,1)+1),[xs(j,end),ys(j,end)],'s','size',2,'color','red');
    end
    if size(xs,1)>1
    strom_tip = insertShape(strom_tip,'Line',reshape([xs(:,end)';ys(:,end)'],1,[]),'LineWidth',1,'color',[150,150,0]);
    end
end
img=[];    
for i=1:60
    img1 = im2double(imread(['RawSplitChannels/',sample,'.tif'],(i-1)*3+2));
    img1 = img1-min(img1(:));
    img1 = img1/max(img1(:));
    img1 = repmat(uint8(img1*255),1,1,3);
    
    img2 = im2double(imread(['RawSplitChannels/segmented_old/chloroplast_',sample,'.tif'],i));
    img2 = uint8(255*img2);
    img2 = cat(3,blank,blank,img2);
    
    img3 = strom_img(:,:,:,i);
    
    img = max(max(max(img1,img2),img3),strom_tip);
    if i==1
        imwrite(img,[sample,'_seg_track_overlay.tiff']);
    else
        imwrite(img,[sample,'_seg_track_overlay.tiff'],'WriteMode','Append');
    end
end

%gifwrite(img,[sample,'_seg_track_overlay.gif']);
%imwrite(strom_tip,[sample,'_tip_motion.png']);