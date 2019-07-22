clear;
close all;
clc;
data_dir = '/data/stromules/output_timeseries/';
lst = dir(data_dir);
cmap1 = hsv(256); %colors for plotting

for fno=3:numel(lst)
    num_slices = numel(dir([data_dir,lst(fno).name]))-2;
    filename = lst(fno).name;
		save_filename = [filename(find(~isspace(filename))),'_new_new'];
		mkdir(save_filename);
    for imno=0:num_slices-1
        tic;
        I = uint8(imread([data_dir,filename,'/t',num2str(imno),'.png'])); %read stromule mask
        I_s = bwmorph(bwareaopen(im2bw(I),'skel',Inf),4); % skeletonize and wipe out small fragments
        RI_s = repmat(I,1,1,3); %make 3 channel for visualizing results
        disp(num2str(imno))
        if imno == 0 %initialize snakes on first frame
            all_stromules = get_snake_pts(I_s);
            all_stromules = evolve_snakes(I_s,all_stromules,5);
        else
            all_stromules = evolve_snakes(I_s,all_stromules,3);
            new_stromules = get_snake_pts(get_remaining(I_s,all_stromules));
            all_stromules = [all_stromules, evolve_snakes(I_s,new_stromules,5)];
        end 
        
				for i=1:numel(all_stromules)
						if ~isempty(all_stromules{i})
								dlmwrite([save_filename,'/',num2str(i),'.txt'],[all_stromules{i}, repmat(imno,size(all_stromules{i},1),1)]','-append');
						end
				end
        % for i=1:numel(all_stromules)
            % npts = all_stromules{i};
            % if ~isempty(npts)
                % %RI_s = insertMarker(RI_s,npts,'size',2);
                % RI_s = insertShape(RI_s,'Line',reshape(npts',1,[]),'LineWidth',3,'color',uint8(255*cmap1(mod(i,256)+1,:)));
            % end
        % end
        % [RI_s,map] = rgb2ind(RI_s,256);
        % if imno == 0
            % imwrite(RI_s,map,[save_filename,'.gif'],'gif', 'Loopcount',inf);
        % else
            % imwrite(RI_s,map,[save_filename,'.gif'],'gif','WriteMode','append');
        % end
        toc;
    end
end