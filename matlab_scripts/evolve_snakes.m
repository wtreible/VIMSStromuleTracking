function all_stromules = evolve_snakes(I_s,all_stromules,num_iter)

parfor i=1:numel(all_stromules)
    npts = all_stromules{i};
    %ppts = npts;
    e1=0;
    if ~isempty(npts)
        for iter=1:num_iter
            %[snake_pnts,e1,e2] = snake(pnts, alpha, beta, max_delta_y, resol_y, max_delta_x, resol_x, feat_img)
            %   pnts          Starting contour. Each row is a [x,y] coordinate.
            %   alpha         Energy contributed by the distance between control points.
            %                 Set to zero if length of slope doesn't matter.
            %   beta          Energy contributed by the curvature of the snake.  Larger
            %                 values of beta cause bends in the snake to have a high cost
            %                 and lead to smoother snakes.
            %   max_delta_y   Max number of pixels to move each contour point vertically
            %   resol_y       Contour points will be moved by multiples of resol_y
            %   max_delta_x   Analog to max_delta_y
            %   resol_x       Analog to resol_y
            %   feat_img      2D-Array of the feature responses in the image.  For example
            %                 it can contain the magnitude of the image gradients
            [npts,e1,e2] = snake(npts,-0.1,0.1,0.001,5,1,5,1,I_s);
        end
        if e2 < 5
            all_stromules{i} = npts;
        else
            all_stromules{i} = [];
        end
    end
end