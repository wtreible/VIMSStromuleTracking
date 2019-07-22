clear;
close all;
clc;
stromule_dir = '/data/inprogress/kolagund/stromule_track/';
chloroplast_dir = '/data/inprogress/kolagund/stromule_track/RawSplitChannels/segmented_old/';
lst = dir([stromule_dir,'*_new']);
dirs = lst([lst.isdir]);

all_changes = [];
all_dl = [];
all_deformation_feat = [];
all_l = [];

for dno=1:numel(dirs)
    sample = dirs(dno).name;
    tracks = dir([stromule_dir, sample, '/*.txt']);
    st_no=0;
    mkdir([stromule_dir, sample, '/updated/']);
    delete([stromule_dir, sample, '/updated/*.txt']);
    for tno=1:numel(tracks)
        track = tracks(tno).name;
        st = dlmread([stromule_dir, sample, '/', track]);
        xs = round(st(1:3:end,:));
        ys = round(st(2:3:end,:));
        zs = round(st(3:3:end,:));
        img_name = ['chloroplast_', sample(1:end-4), '.tif'];
        inds=[0,0];
        max_min=0;
        found=false;
        for i=1:size(zs,1)
            img = imread([chloroplast_dir, img_name],zs(i,1)+1);
            D = bwdist(img);
            vals = D(sub2ind(size(D),[ys(i,1),ys(i,end)], [xs(i,1),xs(i,end)]));
            
            [min_val,ind] = min(vals);
            if any(vals~=min_val)
                found=true;
                inds(ind) = inds(ind)+1;
                if min_val > max_min
                    max_min = min_val;
                end
            end
            %imshow(img);
            %hold on;
            %plot(xs(i,:),ys(i,:),'LineWidth',2);
            %pause(.1);
        end
        if inds(2)>inds(1)
            xs = xs(:,end:-1:1);
            ys = ys(:,end:-1:1);
        end
        xd = xs(:,2:end) - xs(:,1:end-1);
        yd = ys(:,2:end) - ys(:,1:end-1);
        l = sum(sqrt(xd.^2 + yd.^2),2);
        if found && max(l)>=5 && max_min < 15
            disp(st_no);
            dlmwrite([stromule_dir, sample, '/updated/', num2str(st_no),'.txt'],[xs,ys,zs]);
            st_no = st_no+1;
        end
    end
end