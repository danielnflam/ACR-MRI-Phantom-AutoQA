pkg load dicom
pkg load image
pkg load matgeom

clear all
close all

# Settings
centers = [125, 127]; # [row idx, column idx] center pixel of the phantom
sequence = "T2"
switch sequence
  case "T1"
    DICOMDIR_PATH = "Z:\\MRI\\Skyra\\20230706\\T1" # T1
  case "T2"
    DICOMDIR_PATH = "Z:\\MRI\\Skyra\\20230706\\T2" # T2
  case "T2PD"
    DICOMDIR_PATH = "Z:\\MRI\\Skyra\\20230706\\T2" # T2 and PD
endswitch



## DON'T CHANGE BELOW
PixelSpacing = [0.9766, 0.9766]; # [mm]
# Loading Images
[d3d, d3dmeta] = dicomdir_read(DICOMDIR_PATH);
switch sequence
  case "T1"
    d3d = d3d.L;
    d3dmeta = d3dmeta.L;
  case "T2"
    d3d = d3d.L;
    d3dmeta = d3dmeta.L;
  case "T2PD"
    d3d = d3d.H;
    d3dmeta = d3dmeta.H;
endswitch


image = Image3D(d3d, d3dmeta.sliceLocations,...
          d3dmeta.dcmMetadata.PixelSpacing,...
          d3dmeta.dcmMetadata.ImagePositionPatient);

# Action
phantom = ACR_Phantom_QA(centers, image.pixelSpacing)
geom_accuracy_slice1 = phantom.ACR_test_GeometricAccuracy(image.array(:,:,1), false);
saveas(gcf, cat(2, sequence, " geometric accuracy slice 1"), "png")

geom_accuracy_slice5 = phantom.ACR_test_GeometricAccuracy(image.array(:,:,5), true);
saveas(gcf, cat(2, sequence, " geometric accuracy slice 5"), "png")
sprintf("ACR Spec: 190 +/- 2 mm \n")

sideX = 60; sideY = 10; # [pixels] the ROI of the Slice Thickness Accuracy test region
slice_thickness = phantom.ACR_test_SliceThicknessAccuracy(image.array(:,:,1), sideX, sideY);
sprintf("Slice Thickness: %f mm vs ACR: 5+/-0.7 mm\n", slice_thickness)
saveas(gcf, cat(2, sequence, " slice thickness accuracy"), "png")


slice_position_accuracy = phantom.ACR_test_SlicePositionAccuracy(image.array(:,:,1));
sprintf("Slice Position Accuracy Slice 1: %f mm vs ACR: -5<spa<5 mm\n", slice_position_accuracy)
saveas(gcf, cat(2, sequence, " positioning accuracy slice1"), "png")
slice_position_accuracy = phantom.ACR_test_SlicePositionAccuracy(image.array(:,:,11));
sprintf("Slice Position Accuracy Slice 11: %f mm vs ACR: -5<spa<5 mm\n", slice_position_accuracy)
saveas(gcf, cat(2, sequence, " positioning accuracy slice11"), "png")

[slice_uniformity, mean_pixel_intensity] = phantom.ACR_test_SliceUniformity(image.array(:,:,7));
sprintf("Percentage Integral Uniformity: %f %% vs ACR: >87.5 %% for 1.5T, >82 %% for 3T \n", slice_uniformity*100)
saveas (gcf, cat(2, sequence, " uniformity"), "png")

ghosting_ratio = phantom.ACR_test_PercentSignalGhosting(image.array(:,:,7));
sprintf("Ghosting Ratio: %f vs ACR: <0.025 ", ghosting_ratio)
if ghosting_ratio < 0.025; sprintf("PASS \n"); else sprintf("FAIL \n") endif
saveas(gcf, cat(2, sequence, " ghosting ratio"), "png")
