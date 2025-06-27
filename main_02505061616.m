close; clear; clc;

filter = {'Cy5/';'Ftc/'};

for h=1:length(filter)
    f = filter{h};
    directory = dir(append(f,'*.tif'));
    
    condition_bead_list = [];
    condition_radius_list = [];

for i=1:length(directory)
    % get image
    filename = append(f, directory(i).name);
    im = im2double(imread(filename)); 
    
    % create mask
    obj = strel("disk", 3, 8);
    obj2 = strel("disk", 1, 8);
    
    bin = imdilate(imerode(imbinarize(im, "adaptive", 'Sensitivity', 0.57), obj), obj);
    mask = imdilate(bin, obj2);
    
    % Detect circles
    [centers, radii] = imfindcircles(mask, [20 100], 'Sensitivity', 0.8);
    circMask = circles2mask(centers, radii, size(im));
    
    % Remove partial circles and select remaining circles
    cleanMask =  imclearborder(circMask);
    [cleanCenters, cleanRadii] = imfindcircles(cleanMask, [15 100], 'Sensitivity', 0.8);


    % deduce background
    bg = im.*im2double(1-cleanMask);
    bg(bg==0) = NaN;
    
    % interpolate BG and subtract from original image
    bginterp = impaint_nans(bg, 4);
    treated = im - bginterp;
    
    % save treated (Trt) & background (Bkg) images (Pic in 2 folders)
    if i==1
        mkdir(append(f,'Trt'));
        mkdir(append(f,'Bkg'));
    end 
    imwrite(im2uint16(treated), append(f,'Trt/',filename(5:end-4),'_Trt.tif'))
    imwrite(im2uint16(bginterp), append(f,'Bkg/',filename(5:end-4),'_Bkg.tif'))

    Nbeads = length(cleanRadii);
    
    intensities = zeros(Nbeads,1);

    for b=1:Nbeads
        beadMask = double(circles2mask(cleanCenters(b,:), cleanRadii(b), size(im)));
        intensities(b) = sum((beadMask) .* (treated), 'all');
    end
    
    condition_bead_list = [condition_bead_list; intensities];
    condition_radius_list = [condition_radius_list; cleanRadii];
end

BeadNumber = length(condition_bead_list);

figs = gcf;
hi = histogram(condition_bead_list,'BinWidth', 100, 'BinLimits', [0,1800]);
hold on
xlabel('Integrated intensity')
ylabel('Count')
set(gca, "fontsize", 16)
box off
title(append(f,string(datetime('now')),', {\itn_{Beads}} = ',string(BeadNumber)))
xline(median(condition_bead_list),'r','LineWidth',2)
text(median(condition_bead_list)*1.02,max(hi.Values)*-0.02, string(round(median(condition_bead_list))),'color','r')
exportgraphics(figs, append(f(1:end-1), "_comb_int_hist.png"))
close

figs2 = gcf;
hi2 = histogram(condition_radius_list,'BinWidth', 1, 'BinLimits', [20,50]);
hold on
xlabel('Bead radius [px]')
ylabel('Count')
set(gca, "fontsize", 16)
box off
title(append(f,string(datetime('now')),', {\itn_{Beads}} = ',string(BeadNumber)))
xline(median(condition_radius_list),'r','LineWidth',2)
text(median(condition_radius_list)*1.02,max(hi2.Values)*-0.02, string(round(median(condition_radius_list))),'color','r')
exportgraphics(figs2, append(f(1:end-1), "_comb_rad_hist.png"))
close

if ~isempty(condition_bead_list)
A = array2table([condition_radius_list, condition_bead_list]);
A.Properties.VariableNames(1:2) = {'Radius [px]','Integrated intensity [RFU]'};
writetable(A,append(f(1:end-1), "_rad_int.csv"))
end

end