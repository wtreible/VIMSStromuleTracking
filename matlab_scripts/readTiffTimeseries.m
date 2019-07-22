function [ slices ] = readTiffTimeseries(filename)
    fname = filename;
    info = imfinfo(fname);
    num_images = numel(info);
    slice_one = imread(fname, 1, 'Info', info);
    slices = zeros(num_images,size(slice_one,1), size(slice_one, 2),3);
    slices(1,:,:,:) = slice_one;
    for k = 2:num_images
        slices(k,:,:,:) = imread(fname, k, 'Info', info);
    end
end

