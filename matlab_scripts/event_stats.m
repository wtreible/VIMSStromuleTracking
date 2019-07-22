%%%RUN get_stats BEFORE RUNNING THIS SCRIPT.
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
writetable(ET,'EventStats.xls');