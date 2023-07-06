classdef ACR_Phantom_QA
  # Holds functions related to ACR Phantom QA
  properties
    centerX # column idx of phantom centre
    centerY # row idx of phantom centre
    PixelSpacing # [column, row]
  endproperties

  methods
    # NOT using imref
    # TEST FUNCTIONS
    function phantom = ACR_Phantom_QA(centers, PixelSpacing)
      # Initialize with information
      phantom.centerX = centers(2);
      phantom.centerY = centers(1);
      phantom.PixelSpacing = PixelSpacing; # mm
    endfunction

    function geomAccuracy = ACR_test_GeometricAccuracy_localizer(phantom, d3d_slice)
      # ONLY USED FOR LOCALIZER IMAGE_PATH

      # Rough segmentation of water
      [slice_th, threshold] = ACR_Phantom_QA.findPhantomBody(d3d_slice);
      slice_th = imfill(slice_th,"holes");

      # Threshold again; level should be 1/2 mean value of water signal
      water_signals = nonzeros(slice_th.*d3d_slice);
      threshold = mean(water_signals)/2;
      slice_th = d3d_slice > threshold;
      # no hole fill in order to help align shape better

      # Convert to [r,c] data
      [row,col] = find(slice_th);
      d3d_transformed = pca_transform( double([row,col]) , double(d3d_slice) ) ;
      slice_th_transformed = d3d_transformed > threshold;
      slice_th_transformed = imfill(slice_th_transformed, "holes");

      # At center row, find the sagittal length of the phantom
      [row,col] = find(slice_th);
      range_rows = [ min(row), max(row)];
      saggital_idx = []; idxs = [];
      for i = min(range_rows):max(range_rows)
        # find the leftmost and rightmost index at row "i"
        row_i = slice_th_transformed(i,:);
        idx_left = find(row_i,1,"first") ;
        idx_right = find(row_i,1,"last") ;
        saggital_idx = [saggital_idx; idx_right - idx_left];
        idxs = [ idxs ; [idx_left, idx_right] ] ;
      endfor

      geomAccuracy = median(saggital_idx).*phantom.PixelSpacing(1);
      ACR_Phantom_QA.visualize(pca_transform( double([row,col]) , d3d_slice )); hold on;
      plot([median(idxs(:,1)); median(idxs(:,2))] , [120, 120], '-r')
      hold off;
    endfunction

    function geomAccuracy = ACR_test_GeometricAccuracy(phantom, d3d_slice, all_directions)
      # Returns geometric accuracy in milimetres
      # Test geometric accuracy in slice 1
      # Find diameter of sphere in vertical, horizontal & 2 diagonals

      # Rough segmentation of water
      [slice_th, threshold] = ACR_Phantom_QA.findPhantomBody(d3d_slice);
      slice_th = imfill(slice_th,"holes");

      # Level should be 1/2 mean value of water signal
      water_signals = nonzeros(slice_th.*d3d_slice);
      threshold = mean(water_signals)/2;
      slice_th = d3d_slice > threshold;
      # Fill holes for diagonal estimation
      if all_directions
        slice_th = imfill(slice_th,"holes");
      endif

      # Vertical & Horizontal
      centers = [phantom.centerY, phantom.centerX];
      idxs_vert = ACR_Phantom_QA.findCardinalEdgesPhantom(slice_th, "vertical", centers);
      idxs_horz = ACR_Phantom_QA.findCardinalEdgesPhantom(slice_th, "horizontal", centers);
      horz_distance = max(idxs_horz)-min(idxs_horz);
      vert_distance = max(idxs_vert)-min(idxs_vert);
      # Diagonal
      if all_directions
        idx_diags = ACR_Phantom_QA.findDiagonalEdgesPhantom(slice_th, centers);
        diag1_distance = norm(idx_diags(2,:)-idx_diags(1,:));
        diag2_distance = norm(idx_diags(4,:)-idx_diags(3,:));
      endif
      # Calc
      # add 1 pixel spacing for horz & vert, sqrt(2) pixel spacing for diags
      if all_directions
        geomAccuracy = [horz_distance, vert_distance, diag1_distance, diag2_distance]*phantom.PixelSpacing(1) + [1, 1, sqrt(2), sqrt(2)]*phantom.PixelSpacing(1); # mm
      else
        geomAccuracy = [horz_distance, vert_distance]*phantom.PixelSpacing(1) + [1, 1]*phantom.PixelSpacing(1); # mm
      endif

      # Error Reporting:
      disp(" ======= ")
      disp("Geometric Accuracy")
      sprintf("Horizontal: %f  Vertical: %f \n", geomAccuracy(1), geomAccuracy(2))
      if all_directions
        sprintf("DiagonalRALP: %f  DiagonalRPLA: %f \n", geomAccuracy(3), geomAccuracy(4))
      endif
      disp(" ======= ")

      # Visualisation
      ACR_Phantom_QA.visualize(d3d_slice) ; hold on;
      plot([centers(2), centers(2)], [min(idxs_vert); max(idxs_vert)], '-r')
      plot([min(idxs_horz); max(idxs_horz)],[centers(1), centers(1)], '-r')
      if all_directions
        plot([idx_diags(1,1), idx_diags(2,1)], [idx_diags(1,2), idx_diags(2,2)], '-r')
        plot([idx_diags(3,1), idx_diags(4,1)], [idx_diags(3,2), idx_diags(4,2)], '-r')
      endif

    endfunction

    function slice_thickness = ACR_test_SliceThicknessAccuracy(phantom, d3d_slice, sideX, sideY)
      # Slice Thickness Accuracy --> examine the 2 ramps @ the middle of the phantom.
      # 5.0 mm +/- 0.7 mm for ACR.
      [slice_th, threshold] = ACR_Phantom_QA.findPhantomBody(d3d_slice);
      masks = imcomplement(slice_th);
      cc = bwconncomp(masks);
      # find the mask with location (center(2),center(1)) in (row,col)
      for i_mask = 1:length(cc.PixelIdxList)
        [s1,s2] = ind2sub(size(masks), cc.PixelIdxList{i_mask});
        if sum(s1 == phantom.centerY & s2 == phantom.centerX) > 0
          break
        endif
      end
      [s1, s2] = ind2sub(size(masks),cc.PixelIdxList{i_mask});
      mask = false(size(masks)); mask(s1,s2)=true;

      # Crop out the ROI
      wholeROI = imerode(mask, strel("square",5)).*d3d_slice;
      maxInts = nonzeros(max(wholeROI,[],2));
      average_ramp_signal = mean(maxInts);
      # set "display window"
      ROI_thresh_original = wholeROI > average_ramp_signal/2;



      # crop up to 60 pix on each side
      crop_mask = ACR_Phantom_QA.crop_img(mask, phantom.centerX, phantom.centerY, sideX, sideY);
      ROI_thresh = ROI_thresh_original.*crop_mask;

      # Select median option for ramp width
      # BOTTOM
      scores = []; idxs = [];
      for row = phantom.centerY:1:size(ROI_thresh,1)
        idxs_last = find(ROI_thresh(row,:),1,'last');
        idxs_first = find(ROI_thresh(row,:),1,'first');
        diff_idxs = idxs_last - idxs_first;
        if (~isempty(diff_idxs))
          scores = cat(1, scores, diff_idxs);
          idxs = [idxs;[idxs_first, idxs_last]];
        endif
      endfor
      scores = nonzeros(scores);
      score_bottom = median(scores);

      idxs_bottom = [median(idxs(:,1)), median(idxs(:,2))] ; #mean(idxs,1);
      # TOP
      scores = []; idxs = [];
      for row = phantom.centerY-1:-1:1
        idxs_last = find(ROI_thresh(row,:),1,'last');
        idxs_first = find(ROI_thresh(row,:),1,'first');
        diff_idxs = idxs_last - idxs_first;
        if (~isempty(diff_idxs))
          scores = cat(1, scores, diff_idxs);
          idxs = [idxs;[idxs_first, idxs_last]];
        endif
      endfor
      scores = nonzeros(scores);
      score_top = median(scores);
      idxs_top = [median(idxs(:,1)), median(idxs(:,2))] ; # mean(idxs,1);

      # Visualize
      ACR_Phantom_QA.visualize(ROI_thresh_original); hold on;
      plot(idxs_top, [phantom.centerY-3, phantom.centerY-3],'-r')
      plot(idxs_bottom, [phantom.centerY+3, phantom.centerY+3],'-r')



      slice_thickness = phantom.PixelSpacing(1)*0.2*(score_bottom*score_top)/(score_bottom+score_top);
    endfunction
    function slice_position_accuracy = ACR_test_SlicePositionAccuracy(phantom, d3d_slice)
      # Find the slice positioning accuracy by examining the anterior-side crossed wedges
      # Lower slice_position_accuracy is better. < 5 mm for ACR.
      [slice_th,~] = ACR_Phantom_QA.findPhantomBody(d3d_slice);

      idx_upwards = 90; # [px] distance from phantom centre to ~middle of ROI section
      idx_downwards = 10; # [px] distance above central bar (that used for slice thickness accuracy)
      ROI_half_width = 7; #[px]
      ROI_height = 40; #[pix]

      ROI_mask = false(size(d3d_slice));
      ROI_mask( phantom.centerY-idx_upwards:1:phantom.centerY-idx_downwards,...
                 phantom.centerX-ROI_half_width:1:phantom.centerX+ROI_half_width) = true;
      ROI_cut = ROI_mask.*slice_th;
      #imagesc(ROI_cut);

      # From lowest row idx upwards
      lengths_left = []; idxs_left = [];
      for column_idx = phantom.centerX-ROI_half_width:1:phantom.centerX
        idx_left = [find(ROI_cut(:,column_idx), 1, "first"), find(ROI_cut(:,column_idx), 1, "last")];
        diff_idxs = [idx_left(2)-idx_left(1)];
        lengths_left = cat(2, lengths_left, diff_idxs*phantom.PixelSpacing(2));
        idxs_left = [idxs_left;idx_left];
      endfor
      lengths_right = []; idxs_right = [];
      for column_idx = phantom.centerX:1:phantom.centerX+ROI_half_width
        idx_right = [find(ROI_cut(:,column_idx), 1, "first"), find(ROI_cut(:,column_idx), 1, "last")];
        diff_idxs = [idx_right(2) - idx_right(1)];
        lengths_right = cat(2, lengths_right, diff_idxs*phantom.PixelSpacing(2));
        idxs_right = [idxs_right;idx_right];
      endfor
      # Slice position accuracy
      slice_position_accuracy = median(lengths_left) - median(lengths_right); # [mm]

      # Visualization
      ACR_Phantom_QA.visualize(d3d_slice); hold on;
      plot( [phantom.centerX, phantom.centerX] , [median(idxs_left(:,1)), median(idxs_right(:,1))] , '-r')

    endfunction
    function [slice_uniformity, mean_pixel_intensity] = ACR_test_SliceUniformity(phantom, d3d_slice)
      # Find the slice uniformity in slice location 7

      centers = [phantom.centerY, phantom.centerX]; # center of the phantom  (row, col)
      # define a circle of ~200 cm^2 radius
      radius_mm = sqrt(200/pi)*10;
      radius_pixels = floor(radius_mm / phantom.PixelSpacing(1));
      [xx,yy] = meshgrid([1:size(d3d_slice,1)],[1:size(d3d_slice,2)]);
      ROI_mask = (yy - centers(1)).^2 ...
        + (xx - centers(2)).^2 <= radius_pixels.^2;

      # Cut out the ROI
      ROI = ROI_mask.*d3d_slice;
      mean_pixel_intensity = mean(nonzeros(ROI)); # for percentage ghosting in another test

      # Uniformity Testing
      pixel_area = phantom.PixelSpacing(1)*phantom.PixelSpacing(2)/100; #[cm^2]
      [low_ROI] = ACR_Phantom_QA.runUniformityTest(d3d_slice, ROI_mask, centers, radius_pixels, phantom.PixelSpacing, "min");
      [high_ROI] = ACR_Phantom_QA.runUniformityTest(d3d_slice, ROI_mask, centers, radius_pixels, phantom.PixelSpacing, "max");
      slice_uniformity = (1- (high_ROI.mean_signal - low_ROI.mean_signal)/(high_ROI.mean_signal + low_ROI.mean_signal));

      # Visualization
      ACR_Phantom_QA.visualize(d3d_slice);
      hold on; drawEllipse([ phantom.centerX, phantom.centerY, radius_pixels, radius_pixels, 0], 'linewidth',2,'color','red')
      hold on; drawEllipse([ low_ROI.centroid(2), low_ROI.centroid(1), low_ROI.radius_pix, low_ROI.radius_pix, 0], 'linewidth',2,'color','red')
      hold on; drawEllipse([ high_ROI.centroid(2), high_ROI.centroid(1), high_ROI.radius_pix, high_ROI.radius_pix, 0], 'linewidth',2,'color','red')

    endfunction

    function [ghosting_ratio] = ACR_test_PercentSignalGhosting(phantom,d3d_slice)
      # Find percentage signal ghosting in slice location 7

      # Find Phantom Body
      [slice_th, threshold] = ACR_Phantom_QA.findPhantomBody(d3d_slice);
      slice_th = imfill(slice_th,"holes");
      # Level should be 1/2 mean value of water signal
      water_signals = nonzeros(slice_th.*d3d_slice);
      threshold = mean(water_signals)/2;


      # Central ROI Signal
      centers = [phantom.centerY, phantom.centerX]; # center of the phantom  (row, col)
      radius_mm = sqrt(200/pi)*10; # define a circle of ~200 cm^2 radius
      radius_pixels = floor(radius_mm / phantom.PixelSpacing(1));
      [xx,yy] = meshgrid([1:size(d3d_slice,1)],[1:size(d3d_slice,2)]);
      ROI_mask = (yy - centers(1)).^2 ...
        + (xx - centers(2)).^2 <= radius_pixels.^2;
      ROI = ROI_mask.*d3d_slice;
      central_signal = mean(nonzeros(ROI)); # for percentage ghosting in another test

      # Create ellipse ROIs
      ROI_mask = d3d_slice > threshold;
      [top] = ACR_Phantom_QA.runGhostingEllipseTest(d3d_slice, ROI_mask, centers, phantom.PixelSpacing, "top");
      [bottom] = ACR_Phantom_QA.runGhostingEllipseTest(d3d_slice, ROI_mask, centers, phantom.PixelSpacing, "bottom");
      [left] = ACR_Phantom_QA.runGhostingEllipseTest(d3d_slice, ROI_mask, centers, phantom.PixelSpacing, "left");
      [right] = ACR_Phantom_QA.runGhostingEllipseTest(d3d_slice, ROI_mask, centers, phantom.PixelSpacing, "right");

      # Calculations
      ghosting_ratio = abs(((top.signal+bottom.signal)-(left.signal+right.signal))/(2*central_signal));

      # Visualisation
      ACR_Phantom_QA.visualize(d3d_slice);
      hold on; drawEllipse([ phantom.centerX, phantom.centerY, radius_pixels, radius_pixels, 0], 'linewidth',2,'color','red')
      hold on; drawEllipse([top.centroid(2), top.centroid(1), top.X_halfaxis, top.Y_halfaxis, 0], 'linewidth',2,'color','red')
      hold on; drawEllipse([bottom.centroid(2), bottom.centroid(1), bottom.X_halfaxis, bottom.Y_halfaxis, 0], 'linewidth',2,'color','red')
      hold on; drawEllipse([left.centroid(2), left.centroid(1), left.X_halfaxis, left.Y_halfaxis, 0], 'linewidth',2,'color','red')
      hold on; drawEllipse([right.centroid(2), right.centroid(1), right.X_halfaxis, right.Y_halfaxis, 0], 'linewidth',2,'color','red')

    endfunction
  endmethods

  # HELPER FUNCTIONS
  methods (Static = true)
##    function overlayMask(mask, color_rgb)
##      #  DOESN'T WORK YET -- ALPHADATA NOT SUPPORTED BY OCTAVE
##      # Overlay MASK on top of image axis handle H.
##      # COLOR: RGB e.g. green = [0 1 0]
##
##      color_mask = cat(3, ones(size(mask)), ones(size(mask)), ones(size(mask)));
##      for i = 1:3
##        color_mask(:,:,i) = color_mask(:,:,i).*color_rgb(i);
##      endfor
##      hold on
##      h = imshow(color_mask);
##      set(h, 'AlphaData', mask)
##    endfunction
    function h = visualize(d3d_slice)
      d3d_slice = rescale(d3d_slice); # rescale to [0,1]
      img_slice = cat(3,d3d_slice,d3d_slice,d3d_slice);
      figure; h = imshow(img_slice); colormap(gray); axis equal;
    endfunction
    function [output_struct] = runGhostingEllipseTest(d3d_slice, ROI_mask, centers, PixelSpacing, switch_key)
      ROI = ROI_mask.*d3d_slice;
      switch(switch_key)
        case {"top"}
          pixels = [1 : find(ROI(:,centers(2)),1,'first')];
        case {"bottom"}
          pixels = [find(ROI(:,centers(2)),1,'last') : size(d3d_slice,1) ];
        case {"left"}
          pixels = [1 : find(ROI(centers(1),:),1,'first') ];
        case {"right"}
          pixels = [find(ROI(centers(1),:),1,'last') : size(d3d_slice,2) ];
        otherwise
          error("Check Switch Key")
      endswitch
      # Craft Ellipse
      [xx,yy] = meshgrid([1:size(d3d_slice,1)],[1:size(d3d_slice,2)]);
      # Ellipse Area: 1000 mm^2 with 4:1 ratio: area = pi*a*b = pi*4b^2 where b is the shorter axis
      semiminor_axis = sqrt(1000/(4*pi))/PixelSpacing(1); #[pix]
      semimajor_axis = 4*semiminor_axis;
      output_struct = struct;
      switch(switch_key)
       case {"top", "bottom"}
         centroid = [mean(pixels) , centers(2)] ; # [row,col]
         test_ROI_mask = ((yy - centroid(1))/semiminor_axis).^2 ...
            + ((xx - centroid(2))/semimajor_axis).^2 <= 1;
         output_struct.X_halfaxis = semimajor_axis;
         output_struct.Y_halfaxis = semiminor_axis;
       case {"left","right"}
         centroid = [centers(1), mean(pixels)] ; # [row,col]
         test_ROI_mask = ((yy - centroid(1))/semimajor_axis).^2 ...
            + ((xx - centroid(2))/semiminor_axis).^2 <= 1;
         output_struct.X_halfaxis = semiminor_axis;
         output_struct.Y_halfaxis = semimajor_axis;
     endswitch
     output_struct.signal = mean(nonzeros(test_ROI_mask.*d3d_slice));
     output_struct.centroid = centroid;

    endfunction

    function [output_struct] = runUniformityTest(d3d_slice, ROI_mask, centers, radius_pixels, PixelSpacing, switch_key)

      ROI = ROI_mask.*d3d_slice;
      pixel_area = PixelSpacing(1)*PixelSpacing(2)/100; #[cm^2]

      switch(switch_key)
        case {"min","MIN","Min"}
          # Iteratively move intensity level up from minimum
          intensity_level = min(nonzeros(ROI));
          ROI_under_test = ROI >= intensity_level;
        case {"max", "MAX", "Max"}
          intensity_level = max(nonzeros(ROI));
          ROI_under_test = ROI <= intensity_level;
        otherwise
          error("Check Switch Key")
      endswitch
      #ROI_area_original = sum(nonzeros(ROI_under_test)).*pixel_area; #[cm^2],  DEBUG

      flag_continue = true;
      while flag_continue
        switch(switch_key)
          case {"min","MIN","Min"}
            # Iteratively move intensity level up from minimum
            intensity_level += 1;
            ROI_under_test = ROI >= intensity_level;
          case {"max", "MAX", "Max"}
            intensity_level -= 1;
            ROI_under_test = ROI <= intensity_level;
            ROI_under_test = ROI_under_test.*ROI_mask;
        endswitch

        ROI_clusters = abs(ROI_mask - ROI_under_test);
        cc_new = ACR_Phantom_QA.findLargestConnectedComponent(ROI_clusters, 8);

        cluster_area = cc_new.numPixels.*pixel_area; #[cm^2]
        if cluster_area >= 1 # [cm^2]
          # find centroid
          [yy,xx] = ind2sub(cc_new.ImageSize, cc_new.PixelIdxList);
          centroid = [mean(yy), mean(xx)]; #[row, col]
          # stop while loop
          flag_continue = false;
          break
        endif
      endwhile
      #hold on; plot(centroid(2),centroid(1),'xr');

      #check edge
      radius_test_ROI = ceil((sqrt(1/pi)*10)/PixelSpacing(1)); #[pix]
      [centroid] = ACR_Phantom_QA.checkIfSmallROIOutsideLargeROI(centers, centroid, radius_pixels, radius_test_ROI);
      #hold on; plot(centroid(2),centroid(1),'or');
      # draw small ROI
      [xx,yy] = meshgrid([1:size(d3d_slice,1)],[1:size(d3d_slice,2)]);
      test_ROI_mask = (yy - centroid(1)).^2 ...
        + (xx - centroid(2)).^2 <= radius_test_ROI.^2;
      # calculate signal inside test ROI
      output_struct.mean_signal = mean(nonzeros(test_ROI_mask.*d3d_slice));
      output_struct.centroid = centroid;
      output_struct.radius_pix = radius_test_ROI;
    endfunction
    function [centroid] = checkIfSmallROIOutsideLargeROI(centers, centroid, radius_pixels, radius_test_ROI)
      # Check if a small circular ROI, defined by centroid & radius_test_ROI,
      #  is within the large circular ROI defined by centers & radius_pixels.
      desired_radius = radius_pixels-radius_test_ROI;
      if (centroid(1) - centers(1)).^2 + (centroid(2) - centers(2)).^2 > (desired_radius).^2
        # push the centroid away from the ROI edge on a vector defined against centerX & centerY
        vector = [centers(1)-centroid(1), centers(2)-centroid(2)];
        current_radius = norm(vector); #[pix]
        vector = vector./current_radius;#[row, col]
        centroid = centroid + (current_radius - desired_radius).*vector;
      endif
    endfunction
    function [cc_new] = findLargestConnectedComponent(bw, conn)
      cc = bwconncomp (bw, conn);
      # find the largest ccs w/ pixel numbers > mean
      numPixels = [];
      for i = 1:cc.NumObjects
        numPixels = [numPixels, numel(cc.PixelIdxList{i})];
      endfor
      # Select the largest 2 objects
      [numPixels ,idx] = sort(numPixels, 'descend');
      pixelIdxList = cc.PixelIdxList(idx);

      # Implement into cc_new
      cc_new = struct;
      cc_new.numPixels = numPixels(1);
      cc_new.PixelIdxList = pixelIdxList{1};
      cc_new.Connectivity = conn;
      cc_new.ImageSize = cc.ImageSize;
    endfunction
    function [slice_th,threshold] = findPhantomBody(d3d_slice)
      # Find main body of phantom
      t = otsuthresh(imhist(d3d_slice));
      threshold = t*max(d3d_slice(:));
      slice_th = d3d_slice > threshold;
    endfunction

    function idx_of_interest = findCardinalEdgesPhantom(slice_th, opt, centers)
      centerY = centers(1);
      centerX = centers(2);

      if strcmp(opt, "vertical")
        line_test = slice_th(:, centerY);
      elseif strcmp(opt, "horizontal")
        line_test = slice_th(centerX,:);
      endif
      for i = 1:length(line_test)
        if line_test(i) > 0
          break
        endif
      endfor
      top_i = i;
      for i = length(line_test):-1:1
        if line_test(i) > 0
          break
        endif
      endfor
      bottom_i = i;
      idx_of_interest = [top_i, bottom_i];
    endfunction

    function idx_of_interest = findDiagonalEdgesPhantom(slice_th, centers)
      centerY = centers(1);
      centerX = centers(2);
      printf("Assume isotropic in-plane pixels \n")

      # Down-right
      test_centerX = centerX;
      test_centerY = centerY;
      while (slice_th(test_centerY, test_centerX ) > 0)
        test_centerY = test_centerY + 1;
        test_centerX = test_centerX + 1;
      endwhile
      idx_downright = [test_centerX-1, test_centerY-1];
      # Up-right
      test_centerX = centerX;
      test_centerY = centerY;
      while (slice_th(test_centerY, test_centerX ) > 0)
        test_centerY = test_centerY - 1;
        test_centerX = test_centerX + 1;
      endwhile
      idx_upright = [test_centerX-1, test_centerY+1];
      # Down-left
      test_centerX = centerX;
      test_centerY = centerY;
      while (slice_th(test_centerY, test_centerX ) > 0)
        test_centerY = test_centerY + 1;
        test_centerX = test_centerX - 1;
      endwhile
      idx_downleft = [test_centerX+1, test_centerY-1];
      # Up-left
      test_centerX = centerX;
      test_centerY = centerY;
      while (slice_th(test_centerY, test_centerX ) > 0)
        test_centerY = test_centerY - 1;
        test_centerX = test_centerX - 1;
      endwhile
      idx_upleft = [test_centerX+1, test_centerY+1];

      idx_of_interest = [ idx_upleft;
                          idx_downright;
                          idx_downleft;
                          idx_upright];
    endfunction
    function crop_mask = crop_img(mask, centerX, centerY, sideX, sideY)
      sX = [centerX-sideX:1:centerX+sideX]; # columns
      sY = [centerY-sideY:1:centerY+sideY]; # rows
      crop_mask = false(size(mask)); crop_mask(sY,sX) = true;
    endfunction

  endmethods

endclassdef
