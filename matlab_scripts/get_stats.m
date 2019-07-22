clear;
close all;
clc;
data_dir = '/data/inprogress/kolagund/stromule_track/';
lst = dir([data_dir,'*_new']);
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
    tracks = dir([data_dir, sample,'/updated/*.txt']);
    
    event_state=0;
    for tno=1:numel(tracks)
        track = tracks(tno).name;
        st = dlmread([data_dir, sample, '/updated/', track]);
        xs = st(:,1:10);%xs = st(1:3:end,:);
        ys = st(:,11:20);%ys = st(2:3:end,:);
        zs = st(:,21:end);
        xd = xs(:,2:end) - xs(:,1:end-1);
        yd = ys(:,2:end) - ys(:,1:end-1);
        l = sum(sqrt(xd.^2 + yd.^2),2);
        if size(xs,1)>2
            e_s = stcount+1;
            for stno = 1:size(st,1)
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
disp(tot_events);

histogram(cat(1,base_ext_ret{:}),9);
figure; histogram(cat(1,ext_ret{:}),9);
base_ext_ret = cat(1,base_ext_ret{:});
%base_ext_ret(base_ext_ret>-0.15 & base_ext_ret<0.15)=0;
%base_ext_ret = sign(base_ext_ret);
ext_ret = cat(1,ext_ret{:});
%ext_ret(ext_ret>-0.3 & ext_ret<0.3)=0;
%ext_ret = sign(ext_ret);
exp_cont = cat(1,exp_cont{:});

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
writetable(T,'StromuleStats.xls');

% xs = (xs-repmat(mean(xs,2),1,size(xs,2)))./repmat(l(1),size(xs,1),size(xs,2));
% ys = (ys-repmat(mean(ys,2),1,size(ys,2)))./repmat(l(1),size(ys,1),size(ys,2));

% xm = mean(xs,2);
% ym = mean(ys,2);

% if (numel(l)>1) && (all(l>0))
% all_deformation_feat = [all_deformation_feat; [xs(2:end,:) - xs(1:end-1,:), ys(2:end,:) - ys(1:end-1,:)]];
% dl = l(2:end)-l(1:end-1);
% all_changes = [all_changes; sign(dl)];
% all_dl = [all_dl; dl];
% else
% dl=0;
% end
% %disp(dl)
% end
% end
% range = min(abs(all_dl)):1:max(abs(all_dl));
% h1 = histogram(abs(all_dl(all_changes<0)),range);
% h1 = h1.Values';
% h2 = histogram(all_dl(all_changes>0),range);
% h2 = h2.Values';
% histogram(all_changes,3); title('Contraction-Expansion Frequency');
% figure; bar(range(2:end)',[h1,h2]); legend('Contraction Velocity', 'Expansion Velocity'); title('Velocities');

% disp(size(all_deformation_feat));
% K=9;
% [idx, mean_deformation] = kmeans(all_deformation_feat,K,'MaxIter',1000);

% x = -0.5:1/9:0.5;
% y = repmat(0,1,10);
% mx = mean_deformation(:,1:10);
% my = mean_deformation(:,11:20);
% mx = mx + repmat(x,size(mx,1),1);
% my = my + repmat(y,size(my,1),1);
% figure;
% subplot(1,2,1);
% leg = cell(1,K);
% cmap1 = hsv(K);
% for i=1:K
% plot(mx(i,:),my(i,:),'-.','color',cmap1(i,:),'LineWidth',3); hold on;
% leg{i} = num2str(i);
% end
% legend(leg);
% plot(x,y,'-*k','LineWidth',1);
% axis equal;
% hold off;
% title('Deformations');
% subplot(1,2,2); histogram(idx,1:K); title('Distribution')