function TrackStromules(data_dir, output_directory, input_ext, stromule_channel, chloroplast_channel)
% TRACKSTROMULES Tracks stromules segmented from the VIMS Lab Unet Tool
%
% Arguments:
% data_dir (string)______________Full path to the input timeseries data 
%                                segmented from the VIMS Lab Unet Tool
% output_directory (string)______Full path to the folder where output files
%                                will be saved
% input_ext (string)_____________File extension of the input timeseries 
%                                images (Ex. .tif)
% stromule_channel (integer)_____Channel of the segmented image which 
%                                contains stromules (Ex. [R,G,B] = [1,2,3])
% chloroplast_channel (integer)__Channel of the segmented image which
%                                contains chloroplasts (Ex. [R,G,B] = [1,2,3])

% Input Params
% data_dir = '';
% output_directory = '';
% input_ext = '.tif';
% stromule_channel = 2;
% chloroplast_channel = 1;
stromule_channel = str2num(stromule_channel);
chloroplast_channel = str2num(chloroplast_channel);
min_size = 5;

% Mode flags
%(Set everything to true for now...)
track_stromules = true;
correct_stromules = true;
generate_imgs = true;
generate_overlays = true;
generate_stats = true;
generate_figures = true;

% Double check we have a slash on the end...
data_dir = [data_dir, '/'];
output_directory = [output_directory, '/'];

% If input directory is wrong, exit...
if ~exist(data_dir, 'dir')
    disp('Input directory does not exist');
    exit
end

% Create Output Directories if they don't exist
if ~exist(output_directory, 'dir')
   mkdir(output_directory)
end
tracks_dir = [output_directory,'/tracks/'];
if ~exist(tracks_dir, 'dir')
   mkdir(tracks_dir)
end
output_imgs_dir = [output_directory,'/imgs/'];
if ~exist(output_imgs_dir, 'dir')
   mkdir(output_imgs_dir)
end
overlays_dir = [output_directory,'/overlays/'];
if ~exist(overlays_dir, 'dir')
   mkdir(overlays_dir)
end
figures_dir = [output_directory,'/figures/'];
if ~exist(figures_dir, 'dir')
   mkdir(figures_dir)
end
stats_dir = [output_directory, '/stats/'];
if ~exist(stats_dir, 'dir')
   mkdir(stats_dir)
end

% If one of the following modes is true, iterate over all image directories 
% in an input directory (directories of frames)
if track_stromules || correct_stromules || generate_imgs
    lst = dir(data_dir);
    for fno=3:numel(lst)
        num_slices = numel(dir([data_dir,lst(fno).name]))-2;
        filename = lst(fno).name;
        save_filename = [filename(find(~isspace(filename)))];
        if ~exist([tracks_dir, save_filename])
           mkdir([tracks_dir, save_filename]); 
        end
        disp(['Image: ' save_filename])

        % Set the image timer
        image_begin = tic;

        %%%%%%%%%%%%%%%%%%%
        % TRACK STROMULES %
        %%%%%%%%%%%%%%%%%%%
        %=================%
        if track_stromules
            % Clear tracks if they exist already...
            delete([tracks_dir, save_filename, '/*.txt'])
            disp('   >> Tracking stromules...')
            % Iterate over all frames in a image directory
            for imno=0:num_slices-1
                disp(['      >> Time Point: ' num2str(imno)])

                % Set the iteration timer
                iter_begin = tic;

                % Read input image and split stromule channel
                input_fname = [data_dir,'/',filename,'/t',num2str(imno),'_segmented',input_ext];
                I = imread(input_fname);
                I_stromule_ch = I(:,:,stromule_channel);

                % Convert image to BW, remove small objects, skeletonize
                BW = im2bw(I_stromule_ch);
                BW_opened = bwareaopen(BW, min_size);
                BW_skel = bwmorph(BW_opened,'skel',Inf);
                I_s = BW_skel; % (DO NOT) Invert image for snake script

                % Begin iterative Snakes algorithm
                if imno == 0 % Initialize snakes on first frame
                    new_stromules = get_snake_pts(I_s);
                    all_stromules = evolve_snakes(BW_opened,new_stromules,10);
                else
                    all_stromules = evolve_snakes(BW_opened,all_stromules,5);
                    new_stromules = get_snake_pts(get_remaining(I_s,all_stromules));
                    all_stromules = [all_stromules, evolve_snakes(BW_opened,new_stromules,10)];
                end 

                % Append snake coordinates to the tracks output text files
                for i=1:numel(all_stromules)
                    curr_stromules = all_stromules{i};
                    if ~isempty(curr_stromules)
                        dlmwrite([tracks_dir, save_filename, '/', num2str(i), '.txt'],[curr_stromules, repmat(imno,size(curr_stromules,1),1)]','-append');
                    end
                end

                % Get iteration time and print
                iter_time = toc(iter_begin);
                disp(['      >> Time Point Duration (s): ' num2str(iter_time)])
            end
        end
        %=====================%
        %%%%%%%%%%%%%%%%%%%%%%%
        % END TRACK STROMULES %
        %%%%%%%%%%%%%%%%%%%%%%%  
        
        %%%%%%%%%%%%%%%%%%%%%
        % CORRECT STROMULES %
        %%%%%%%%%%%%%%%%%%%%%
        %===================%
        if correct_stromules
            % Fix snake coordinates and make corrected track dir
            disp('   >> Fixing stromule endpoints...')
            corrected_dir = [tracks_dir, save_filename, '/corrected/'];
            if ~exist(corrected_dir, 'dir')
                mkdir(corrected_dir);
            end
            delete([tracks_dir, save_filename, '/corrected/*.txt']);

            % Initialize new stromule ids to start at 1
            st_no = 1;

            % For each track for this image...
            tracks = dir([tracks_dir, save_filename, '/*.txt']);
            for tno=1:numel(tracks)
                trackfname = tracks(tno).name;
                disp(['      >> Analyzing Track: ', trackfname])
                st = dlmread([tracks_dir, save_filename, '/', trackfname]);
                xs = round(st(1:3:end,:));
                ys = round(st(2:3:end,:));
                zs = round(st(3:3:end,:));

                % Correct orientation of stromules
                inds=[0,0];
                max_min=0;
                found=false;
                for i=1:size(zs,1)
                    input_fname = [data_dir,'/',filename,'/t',num2str(zs(i,1)),'_segmented',input_ext];
                    I2 = imread(input_fname);
                    I_chloro_ch = I2(:,:,chloroplast_channel);
                    
                    % Determine distances to chloroplasts
                    D = bwdist(I_chloro_ch);
                    vals = D(sub2ind(size(D),[ys(i,1),ys(i,end)], [xs(i,1),xs(i,end)]));
                    
                    % Find the largest minimum
                    [min_val,ind] = min(vals);
                    if any(vals~=min_val)
                        found=true;
                        inds(ind) = inds(ind)+1;
                        if min_val > max_min
                            max_min = min_val;
                        end
                    end
                end
                if inds(2)>inds(1)
                    xs = xs(:,end:-1:1);
                    ys = ys(:,end:-1:1);
                end
                xd = xs(:,2:end) - xs(:,1:end-1);
                yd = ys(:,2:end) - ys(:,1:end-1);
                l = sum(sqrt(xd.^2 + yd.^2),2);
                if found && max(l) >= 5 && max_min < 15
                    disp(['      >> Saving as Corrected Stromule Track: ', num2str(st_no)]);
                    corrected_tracks = zeros(3*size(zs,1), size(zs,2));
                    corrected_tracks(1:3:end,:) = xs;
                    corrected_tracks(2:3:end,:) = ys;
                    corrected_tracks(3:3:end,:) = zs;
                    dlmwrite([corrected_dir, num2str(st_no),'.txt'],corrected_tracks);
                    st_no = st_no+1;
                end
            end
        end
        %=======================%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % END CORRECT STROMULES %
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % GENERATE OUTPUT IMAGES %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %========================%
        if generate_imgs
            % Initialize output images if generate_imgs is true
            first_frame = imread([data_dir,'/',filename,'/t0_segmented',input_ext]);
            blank_size = size(first_frame);
            blank_size = blank_size(1:2);
            blank = uint8(zeros(blank_size));
            strom_img = repmat(blank,1,1,3,num_slices);
            strom_tip = repmat(blank,1,1,3);
            clear first_frame; clear blank_size; clear blank;

            % If corrected tracks exist...
            if exist([output_directory, '/tracks/', save_filename, '/corrected/'], 'dir')
                tracks_dirname = ['/tracks/', save_filename, '/corrected/'];
            else
                tracks_dirname = ['/tracks/', save_filename, '/'];
            end

            % Set tracks directory and iterate over tracks
            tracks = dir([output_directory, tracks_dirname, '*.txt']);
            for tno=1:numel(tracks)
                track = tracks(tno).name;
                st = dlmread([output_directory, tracks_dirname, track]);
                xs = st(1:3:end, :);
                ys = st(2:3:end, :);
                zs = st(3:3:end, :);

                % Loop over frames this track exists and draw the tracks and
                % tips on the output strom_img
                for j = 1:size(zs,1)
                    strom_img(:,:,:,zs(j,1)+1) = insertShape(strom_img(:,:,:,zs(j,1)+1),'Line',reshape([xs(j,:);ys(j,:)],1,[]),'LineWidth',3,'color',[0,255,0]);
                    strom_img(:,:,:,zs(j,1)+1) = insertMarker(strom_img(:,:,:,zs(j,1)+1),[xs(j,1),ys(j,1)],'s','size',2,'color','blue');
                    strom_img(:,:,:,zs(j,1)+1) = insertMarker(strom_img(:,:,:,zs(j,1)+1),[xs(j,end),ys(j,end)],'s','size',2,'color','red');
                end
                if size(xs,1)>1
                    strom_tip = insertShape(strom_tip,'Line',reshape([xs(:,end)';ys(:,end)'],1,[]),'LineWidth',1,'color',[150,150,0]);
                end
            end

            % Write tip image and stromule image
            imwrite(strom_tip, [output_imgs_dir, save_filename, '_tip_motion.tif']);
            for imno=1:num_slices
                out_img = strom_img(:,:,:,imno);
                if imno==1
                    imwrite(out_img,[output_imgs_dir, save_filename, '_stromules.tif']);
                else
                    imwrite(out_img,[output_imgs_dir, save_filename, '_stromules.tif'],'WriteMode','Append');
                end
            end
        end
        %============================%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % END GENERATE OUTPUT IMAGES %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Calculate final runtime
        image_time = toc(image_begin);
        disp(['Image Elapsed Time: ' num2str(image_time)])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET STROMULE STATISTICS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=========================%
if generate_stats
    disp('Writing statistics for dataset...');
    lst = dir(tracks_dir);
    dirs = lst([lst.isdir]);

    all_changes = [];
    all_dl = [];
    all_deformation_feat = [];
    all_l = [];

    filename=[];
    stromuleId=[];
    stromuleXY=[];
    frameNo=[];
    stromuleBaseAngle=[];
    stromuleTipAngle=[];
    stromuleLength=[];
    baseMovementAngle=[];
    tipMovementAngle=[];
    tipSpeed=[];
    tipAcc=[];
    baseSpeed=[];
    baseAcc=[];
    curvature=[];
    ext_ret=[];
    base_ext_ret=[];
    exp_cont=[];
    num_strmo_per_image=[];
    tot_events=0;
    stcount=0;
    event_boundaries=[];
    e_s = 0;
    e_e = 0;
    for dno=1:numel(dirs)
        sample = dirs(dno).name;
        
        if exist([tracks_dir, sample, '/corrected/'], 'dir')
            corrected_dirname = [sample, '/corrected/'];
        else
            corrected_dirname = [sample, '/'];
        end
        %tracks = dir([tracks_dir, sample, '/corrected/*.txt']);
        tracks = dir([tracks_dir, corrected_dirname, '*.txt']);

        event_state=0;
        for tno=1:numel(tracks)
            track = tracks(tno).name;
            %st = dlmread([tracks_dir, sample, '/corrected/', track]);
            st = dlmread([tracks_dir, corrected_dirname, track]);
            
            xs = st(1:3:end, :);
            ys = st(2:3:end, :);
            zs = st(3:3:end, :);

            xd = xs(:,2:end) - xs(:,1:end-1);
            yd = ys(:,2:end) - ys(:,1:end-1);
            l = sum(sqrt(xd.^2 + yd.^2),2);
            
            xs = (xs-repmat(mean(xs,2),1,size(xs,2)))./repmat(l(1),size(xs,1),size(xs,2));
            ys = (ys-repmat(mean(ys,2),1,size(ys,2)))./repmat(l(1),size(ys,1),size(ys,2));

            xm = mean(xs,2);
            ym = mean(ys,2);

            if (numel(l)>1) && (all(l>0))
                all_deformation_feat = [all_deformation_feat; [xs(2:end,:) - xs(1:end-1,:), ys(2:end,:) - ys(1:end-1,:)]];
                dl = l(2:end)-l(1:end-1);
                all_changes = [all_changes; sign(dl)];
                all_dl = [all_dl; dl];
            else
                dl=0;
            end
            
            if size(xs,1)>2
                e_s = stcount+1;
                for stno = 1:size(xs,1)
                    stcount=stcount+1;
                    filename{end+1} = sample;
                    stromuleId{end+1} = ['ST', sprintf( '%06d', tno )];
                    ppx = spline(0:9,xs(stno,:));
                    ppy = spline(0:9,ys(stno,:));
                    ppdx = fnder(ppx);
                    ppdxx = fnder(ppdx);
                    ppdy = fnder(ppy);
                    ppdyy = fnder(ppdy);
                    curv = abs((ppval(ppdx,1:8).*ppval(ppdyy,1:8) - ppval(ppdxx,1:8).*ppval(ppdy,1:8))./((ppval(ppdx,1:8).^2 + ppval(ppdy,1:8).^2).^1.5));
                    curv(isnan(curv))=0;
                    curvature{end+1} = curv;
                    stromuleXY{end+1} = reshape([xs(stno,:);ys(stno,:)],1,[]);
                    frameNo{end+1} = zs(stno,1);
                    baseLine = [xs(stno,3)-xs(stno,1),ys(stno,3)-ys(stno,1)];
                    stromuleBaseAngle{end+1} = baseLine/norm(baseLine);
                    tipLine = [xs(stno,end)-xs(stno,end-2),ys(stno,end)-ys(stno,end-2)];
                    stromuleTipAngle{end+1} = tipLine/norm(tipLine);
                    stromuleLength{end+1} = l(stno);
                    if stno>2
                        baseMovement = [xs(stno,1)-xs(stno-1,1),ys(stno,1)-ys(stno-1,1)];
                        baseMovementAngle{end+1} = baseMovement/norm(baseMovement);
                        tipMovement = [xs(stno,end)-xs(stno-1,end),ys(stno,end)-ys(stno-1,end)];
                        tipMovementAngle{end+1} = tipMovement/norm(tipMovement);
                        tipSpeed{end+1} = norm(tipMovement);
                        tipAcc{end+1} = tipSpeed{end} - tipSpeed{end-1};
                        baseSpeed{end+1} = norm(baseMovement);
                        baseAcc{end+1} = baseSpeed{end} - baseSpeed{end-1};
                        ext_ret{end+1} = stromuleTipAngle{end-1}*tipMovementAngle{end}';
                        base_ext_ret{end+1} = stromuleBaseAngle{end-1}*baseMovementAngle{end}';
                        exp_cont{end+1} = stromuleLength{end} - stromuleLength{end-1};
                    elseif stno>1
                        baseMovement = [xs(stno,1)-xs(stno-1,1),ys(stno,1)-ys(stno-1,1)];
                        baseMovementAngle{end+1} = baseMovement/norm(baseMovement);
                        tipMovement = [xs(stno,end)-xs(stno-1,end),ys(stno,end)-ys(stno-1,end)];
                        tipMovementAngle{end+1} = tipMovement/norm(tipMovement);
                        tipSpeed{end+1} = norm(tipMovement);
                        tipAcc{end+1} = nan;
                        baseSpeed{end+1} = norm(baseMovement);
                        baseAcc{end+1} = nan;                    
                        ext_ret{end+1} = stromuleTipAngle{end-1}*tipMovementAngle{end}';
                        base_ext_ret{end+1} = stromuleBaseAngle{end-1}*baseMovementAngle{end}';
                        exp_cont{end+1} = stromuleLength{end} - stromuleLength{end-1};
                    else
                        baseMovementAngle{end+1} = [nan, nan];
                        tipMovementAngle{end+1} = [nan, nan];
                        tipSpeed{end+1} = nan;
                        tipAcc{end+1} = nan;
                        baseSpeed{end+1} = nan;
                        baseAcc{end+1} = nan;
                        ext_ret{end+1} = nan;
                        base_ext_ret{end+1} = nan;
                        exp_cont{end+1} = nan;
                    end
                    if stno>1
                        chng = ext_ret{end}*ext_ret{end-1};
                        event_det = ((~isnan(chng) & (chng <= 0)) & tipSpeed{end}>0);
                        if event_det
                            e_e = stcount-1;
                            event_boundaries{end+1}=[e_s, e_e];
                            e_s = stcount-1;
                        end
                        tot_events = tot_events+double(event_det);
                    end
                end
                if ~event_det
                   e_e = stcount;
                   event_boundaries{end+1}=[e_s, e_e]; 
                end
            end
        end
       num_strmo_per_image{end+1}=stcount;
    end

    base_ext_ret = cat(1,base_ext_ret{:});
    ext_ret = cat(1,ext_ret{:});
    exp_cont = cat(1,exp_cont{:});
    
    %histogram(base_ext_ret,9);
    %outfig = figure('visible','off'); 
    %histogram(ext_ret,9);
    %saveas(outfig,[figures_dir, 'ExtensionRetraction.png'])
    
    T = table(filename', cat(1,stromuleId{:}), num2str(cat(1,stromuleXY{:})),...
        cat(1,stromuleBaseAngle{:}), cat(1,baseMovementAngle{:}), cat(1,baseSpeed{:}), cat(1,baseAcc{:}),...
        cat(1,stromuleTipAngle{:}), cat(1,tipMovementAngle{:}), cat(1,tipSpeed{:}), cat(1,tipAcc{:}),...
        cat(1,stromuleLength{:}),cat(1,curvature{:}),...
        ext_ret, base_ext_ret, exp_cont,...
        'VariableNames',{'Filename';'StromuleId';'StromuleXY';...
        'StromuleBaseAngle';'StromuleMovementAngle';'StromuleBaseSpeed';'StromuleBaseAcc';...
        'StromuleTipAngle';'StromuleTipMovementAngle';'StromuleTipSpeed';'StromuleTipAcc';...
        'StromuleLength';'StromuleCurvature';...
        'StromuleTipExtRet'; 'StromuleBaseExtRet'; 'StromuleExpCont'});
    writetable(T,[stats_dir, 'StromuleStats.xls']);
    
    %%%%%%%%%%%
    % Figures %
    %%%%%%%%%%%
    if generate_figures
        
        range = min(abs(all_dl)):1:max(abs(all_dl));
        h1 = histogram(abs(all_dl(all_changes<0)),range);
        h1 = h1.Values';
        h2 = histogram(all_dl(all_changes>0),range);
        h2 = h2.Values';
        
        % Figure 1
        outfig = figure('visible','off');
        histogram(all_changes,3); title('Contraction-Expansion Frequency');
        saveas(outfig,[figures_dir, 'ContractionExpansionFrequency.png'])
        
        % Figure 2
        outfig = figure('visible','off'); 
        bar(range(2:end)',[h1,h2]); legend('Contraction Velocity', 'Expansion Velocity'); title('Velocities');
        saveas(outfig,[figures_dir, 'Velocities.png'])

        K=min(size(all_deformation_feat,1),9);
        [idx, mean_deformation] = kmeans(all_deformation_feat,K,'MaxIter',1000);
        x = -0.5:1/9:0.5;
        y = repmat(0,1,10);
        mx = mean_deformation(:,1:10);
        my = mean_deformation(:,11:20);
        mx = mx + repmat(x,size(mx,1),1);
        my = my + repmat(y,size(my,1),1);
        
        % Figure 3
        outfig = figure('visible','off'); 
        subplot(1,2,1);
        leg = cell(1,K);
        cmap1 = hsv(K);
        hold on;
        for i=1:K
            plot(mx(i,:),my(i,:),'-.','color',cmap1(i,:),'LineWidth',3); hold on;
            leg{i} = num2str(i);
        end
        legend(leg);
        plot(x,y,'-*k','LineWidth',1);
        axis equal;
        hold off;
        title('Deformations');
        subplot(1,2,2); histogram(idx,1:K); title('Distribution')
        saveas(outfig,[figures_dir, 'Deformations.png'])
    end
    
    %%%%%%%%%%
    % Events %
    %%%%%%%%%%
    disp('Generating event statistics for dataset...')
    eId=1;
    Name=[];
    StId=[];
    eventType=[];
    eventId=[];
    eventLength=[];
    eventStart=[];
    eventStop=[];
    avgLength=[];
    minLength=[];
    maxLength=[];
    avgTipSpeed=[];
    minTipSpeed=[];
    maxTipSpeed=[];
    avgTipAcc=[];
    minTipAcc=[];
    maxTipAcc=[];
    avgBaseSpeed=[];
    minBaseSpeed=[];
    maxBaseSpeed=[];
    avgBaseAcc=[];
    minBaseAcc=[];
    maxBaseAcc=[];
    avgCurvature=[];
    minCurvature=[];
    maxCurvature=[];
    avgBaseMovementCorelation=[];
    minBaseMovementCorelation=[];
    maxBaseMovementCorelation=[];
    avgTipMovementCorelation=[];
    minTipMovementCorelation=[];
    maxTipMovementCorelation=[];
    avgLengthChange=[];
    minLengthChange=[];
    maxLengthChange=[];

    for i=1:numel(event_boundaries)
        eb = event_boundaries{i};
        if eb(2)-eb(1) > 2
            Name{end+1} = T{eb(1),1}{1};
            StId{end+1} = T{eb(1),2};
            eventId{end+1} = ['E', sprintf( '%06d', eId )]; eId=eId+1;
            eventtyp = T{eb(1)+1:eb(2),14};
            eventtyp = eventtyp(~isnan(eventtyp));
            if ~isempty(eventtyp)
                eventType{end+1} = sign(eventtyp(1));
            else
                eventType{end+1}=nan;
            end
            eventLength{end+1} = eb(2)-eb(1);
            eventStart{end+1} = eb(1);
            eventStop{end+1} = eb(2);

            lengths = T{eb(1):eb(2),12};
            avgLength{end+1} = mean(lengths);
            minLength{end+1} = min(lengths);
            maxLength{end+1} = max(lengths);

            tipSpeeds = T{eb(1)+1:eb(2),10};
            tipSpeeds = tipSpeeds(~isnan(tipSpeeds));
            baseSpeeds = T{eb(1)+1:eb(2),6};
            baseSpeeds = baseSpeeds(~isnan(baseSpeeds));
            if ~isempty(tipSpeeds)
                avgTipSpeed{end+1} = mean(tipSpeeds);
                minTipSpeed{end+1} = min(tipSpeeds);
                maxTipSpeed{end+1} = max(tipSpeeds);
            else
                avgTipSpeed{end+1} = nan;
                minTipSpeed{end+1} = nan;
                maxTipSpeed{end+1} = nan;
            end
            if ~isempty(baseSpeeds)
                avgBaseSpeed{end+1} = mean(baseSpeeds);
                minBaseSpeed{end+1} = min(baseSpeeds);
                maxBaseSpeed{end+1} = max(baseSpeeds);
            else
                avgBaseSpeed{end+1} = nan;
                minBaseSpeed{end+1} = nan;
                maxBaseSpeed{end+1} = nan;
            end

            tipAcc = T{eb(1)+1:eb(2),11};
            tipAcc = tipAcc(~isnan(tipAcc));
            baseAcc = T{eb(1)+1:eb(2),7};
            baseAcc = baseAcc(~isnan(baseAcc));
            if ~isempty(tipAcc)
                avgTipAcc{end+1} = mean(tipAcc);
                minTipAcc{end+1} = min(tipAcc);
                maxTipAcc{end+1} = max(tipAcc);
            else
                avgTipAcc{end+1} = nan;
                minTipAcc{end+1} = nan;
                maxTipAcc{end+1} = nan;
            end
            if ~isempty(baseAcc)
                avgBaseAcc{end+1} = mean(baseAcc);
                minBaseAcc{end+1} = min(baseAcc);
                maxBaseAcc{end+1} = max(baseAcc);
            else
                avgBaseAcc{end+1} = nan;
                minBaseAcc{end+1} = nan;
                maxBaseAcc{end+1} = nan;
            end

            curvature = T{eb(1):eb(2),13};
            avgCurvature{end+1} = mean(curvature);
            minCurvature{end+1} = min(curvature);
            maxCurvature{end+1} = max(curvature);

            basemovement = T{eb(1)+1:eb(2),15};
            basemovement = basemovement(~isnan(basemovement));
            tipmovement = T{eb(1)+1:eb(2),14};
            tipmovement = tipmovement(~isnan(tipmovement));

            if ~isempty(basemovement)
                avgBaseMovementCorelation{end+1} = mean(basemovement);
                minBaseMovementCorelation{end+1} = min(basemovement);
                maxBaseMovementCorelation{end+1} = max(basemovement);
            else
                avgBaseMovementCorelation{end+1} = nan;
                minBaseMovementCorelation{end+1} = nan;
                maxBaseMovementCorelation{end+1} = nan;
            end
            if ~isempty(tipmovement)
                avgTipMovementCorelation{end+1} = mean(tipmovement);
                minTipMovementCorelation{end+1} = min(tipmovement);
                maxTipMovementCorelation{end+1} = max(tipmovement);
            else
                avgTipMovementCorelation{end+1} = nan;
                minTipMovementCorelation{end+1} = nan;
                maxTipMovementCorelation{end+1} = nan;
            end

            lengthchange = T{eb(1)+1:eb(2),16};
            avgLengthChange{end+1} = mean(lengthchange);
            minLengthChange{end+1} = min(lengthchange);
            maxLengthChange{end+1} = max(lengthchange);
        end
    end
    ET = table(Name', cat(1,eventId{:}), cat(1,eventType{:}), cat(1,eventLength{:}), cat(1,StId{:}),...
        cat(1,eventStart{:}), cat(1,eventStop{:}),...
        cat(1,avgLength{:}), cat(1,minLength{:}), cat(1,maxLength{:}),...
        cat(1,avgTipSpeed{:}), cat(1,minTipSpeed{:}), cat(1,maxTipSpeed{:}),...
        cat(1,avgTipAcc{:}), cat(1,minTipAcc{:}), cat(1,maxTipAcc{:}),...
        cat(1,avgBaseSpeed{:}), cat(1,minBaseSpeed{:}), cat(1,maxBaseSpeed{:}),...
        cat(1,avgBaseAcc{:}), cat(1,minBaseAcc{:}), cat(1,maxBaseAcc{:}),...
        cat(1,avgCurvature{:}), cat(1,minCurvature{:}), cat(1,maxCurvature{:}),...
        cat(1,avgTipMovementCorelation{:}), cat(1,minTipMovementCorelation{:}), cat(1,maxTipMovementCorelation{:}),...
        cat(1,avgBaseMovementCorelation{:}), cat(1,minBaseMovementCorelation{:}), cat(1,maxBaseMovementCorelation{:}),...
        cat(1,avgLengthChange{:}), cat(1,minLengthChange{:}), cat(1,maxLengthChange{:}),...
        'VariableNames',{'Filename'; 'EventId'; 'EventType'; 'EventLength'; 'StromuleId';...
        'EventStart'; 'EventStop';...
        'MeanLength';'MinLength';'MaxLength';...
        'MeanTipSpeed';'MinTipSpeed';'MaxTipSpeed';...
        'MeanTipAcc';'MinTipAcc';'MaxTipAcc';...
        'MeanBaseSpeed';'MinBaseSpeed'; 'MaxBaseSpeed';...
        'MeanBaseAcc'; 'MinBaseAcc'; 'MaxBaseAcc';...
        'MeanCurvature'; 'MinCurvature'; 'MaxCurvature';...
        'MeanTipMovementCorelation'; 'MinTipMovementCorelation'; 'MaxTipMovementCorelation';...
        'MeanBaseMovementCorelation'; 'MinBaseMovementCorelation'; 'MaxBaseMovementCorelation';...
        'MeanLengthChange'; 'MinLengthChange'; 'MaxLengthChange'});
    writetable(ET,[stats_dir, 'EventStats.xls']);
end
%=============================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END GET STROMULE STATISTICS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

