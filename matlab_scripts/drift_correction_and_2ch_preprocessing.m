chloroplast_dir = '/data/inprogress/kolagund/stromule_track/RawSplitChannels/';
lst = dir([chloroplast_dir,'*.tif']);
for j=1:3%numel(lst)
    figure;
    for ch=1:2
        save_name = [lst(j).name(1:end-4),'_ch1.tif'];
        img1 = im2double(imread([chloroplast_dir,lst(j).name],ch));
        img1 = (img1-min(img1(:)))/(max(img1(:))-min(img1(:)));
        x=[0];
        y=[0];
        for i=2:60
            img2 = im2double(imread([chloroplast_dir,lst(j).name],(i-1)*3+ch));
            img2 = (img2-min(img2(:)))/(max(img2(:))-min(img2(:)));
            tform = imregcorr(img2,img1,'translation');
            %disp(tform.T);
            img1 = imwarp(img2,tform);
            
            x=[x;x(end)+tform.T(end,1)];
            y=[y;y(end)+tform.T(end,2)];
            
            %imwrite(uint8(img2*255),[chloroplast_dir,'Motion_corrected/',save_name],'WriteMode','Append');
        end
        hold on;
        quiver(x(1:end-1),y(1:end-1),x(2:end)-x(1:end-1),y(2:end)-y(1:end-1));
        plot(x([1,end]),y([1,end]),'*'); axis equal;
    end
end

%% 
data_dir = '/data/inprogress/kolagund/stromules/2CHTraining/masks/T30/8-1-2016/';
lst = dir([data_dir,'*.tiff']);
for i=1:numel(lst)
img = imread([data_dir,lst(i).name]);
img1 = im2bw(img(:,:,1));
img2 = im2bw(img(:,:,2));

D=bwdist(~img1);
D=-D;
D(~img1)=Inf;
L = watershed(D);
L(~img1)=0;
img1 = L>0;

D = bwdist(img2);
[r,c] = find(D<15);
img1 = bwselect(img1,c,r,8);

D = bwdist(img1);
[r,c] = find(D<15);
img2 = bwselect(img2,c,r,8);

img12 = imclose(img1|img2,ones(15));

img2 = img2 | (~imclose(img1,ones(5)) & img12) & ~im2bw(imclose(img(:,:,1),ones(25)));
img1 = imclose(img1,ones(5)) | (~img2 & img12) & ~im2bw(img(:,:,2));

img(:,:,1) = uint8(255*img1);
img(:,:,2) = uint8(255*img2);

imwrite(img,[data_dir,lst(i).name]);

end

%% 
data_dir = '/data/inprogress/kolagund/stromules/Labels/';
lst = dir(data_dir);
for i=3:numel(lst)
    mkdir([data_dir,'masks/',lst(i).name]);
    ims = dir([data_dir,lst(i).name,'/CH2-T1-CP-Body/*.tif']);
    for j=1:numel(ims)
        save_name = strsplit(ims(j).name,'.');
        save_name = [save_name{1},'.tiff'];
        
        img2 = imbinarize(imread([data_dir,lst(i).name,'/CH2-T1/',ims(j).name]));
        img2 = imfill(img2,'holes');
        
        img1 = imbinarize(imread([data_dir,lst(i).name,'/CH2-T1-CP-Body/',ims(j).name]));
        img1 = imfill(img1,'holes');
        
%         D=bwdist(~img1);
%         D=-D;
%         D(~img1)=Inf;
%         L = watershed(D);
%         L(~img1)=0;
%         img1 = L>0;
        
        D = bwdist(img2);
        [r,c] = find(D<5);
        img1 = bwselect(img1,c,r,8);
        
        D = bwdist(img1);
        [r,c] = find(D<5);
        img2 = bwselect(img2,c,r,8);
        
        %img12 = imclose(img1|img2,ones(11));
        
        %img1 = bwareaopen(img12 & ~img2, 100);
        %img2 = img12 & ~img1;
        
        img2 = imdilate(bwmorph(img2,'thin',Inf),ones(3));
        
        img = uint8(255*cat(3,img1,img2,false(size(img1))));
        
        %imshow(img);
        %pause;
        imwrite(img,[data_dir,'masks/',lst(i).name,'/',save_name]);
    end
end
%% 

%% 
data_dir = '/data/inprogress/kolagund/stromules/Labels/';
lst = dir(data_dir);
for i=3:numel(lst)
    mkdir([data_dir,'masks/',lst(i).name]);
    ims = dir([data_dir,lst(i).name,'/CH2-T1-CP-Body/*.tif']);
    for j=1:numel(ims)
        save_name = strsplit(ims(j).name,'.');
        save_name = [save_name{1},'.tiff'];

        img1 = imbinarize(imread([data_dir,lst(i).name,'/CH2-T1-CP-Body/',ims(j).name]));
        img1 = imfill(img1,'holes');
        
        img2 = imbinarize(imread([data_dir,lst(i).name,'/CH2-T1/',ims(j).name]));
        img2 = imfill(img2,'holes');
        
        img2 = imdilate(thin(img2),ones(3));
        
        img = uint8(255*cat(3,img1,img2,false(size(img1))));
        
        %imshow(img); pause;

        imwrite(img,[data_dir,'masks/',lst(i).name,'/',save_name]);
    end
end
%% 

data_dir = '/data/inprogress/kolagund/stromules/';

lst = dir([data_dir,'tch_output_body/*_predicted.png']);

for i=1:numel(lst)
    img = max(imread([data_dir,'tch_output_body/',lst(i).name]), imread([data_dir, 'tch_output_stromule/', lst(i).name]));
    imwrite(img,[data_dir,'tch_combined/',lst(i).name]);
end


