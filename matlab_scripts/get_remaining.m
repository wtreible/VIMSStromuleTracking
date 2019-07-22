function I_s = get_remaining(I_s, all_stromules)
for i=1:numel(all_stromules)
    xy = all_stromules{i};
    if ~isempty(xy)
        I_s = ~bwselect(I_s,xy(:,1),xy(:,2),8) & I_s;
    end
end
end