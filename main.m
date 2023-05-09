pkg load dicom
pkg load image
pkg load matgeom

clear all
close all

# Settings
DICOMDIR_PATH = "Z:\\MRI\\Images\\20230504\\MRT1\\PZPDWI4U\\L3VUFO5S" # T1
#DICOMDIR_PATH = "Z:\\MRI\\Images\\20230504\\MRT2\\PZPDWI4U\\TRWPWDEW" # T2


centers = [122, 128]; # [row idx, column idx] center pixel of the phantom
PixelSpacing = [0.9766, 0.9766]; # [mm]

# Loading Images
[d3d, d3dmeta] = dicomdir_read(DICOMDIR_PATH);


image = Image3D(d3d, d3dmeta.sliceLocations,...
          d3dmeta.dcmMetadata.PixelSpacing,...
          d3dmeta.dcmMetadata.ImagePositionPatient);

# Action
phantom = ACR_Phantom_QA(centers, image.pixelSpacing)
geom_accuracy_slice1 = phantom.ACR_test_GeometricAccuracy(image.array(:,:,1), false);
geom_accuracy_slice5 = phantom.ACR_test_GeometricAccuracy(image.array(:,:,5), true);
sprintf("ACR Spec: 190 +/- 2 mm \n")

sideX = 60; sideY = 10; # [pixels] the ROI of the Slice Thickness Accuracy test region
slice_thickness = phantom.ACR_test_SliceThicknessAccuracy(image.array(:,:,1), sideX, sideY);
sprintf("Slice Thickness: %f mm vs ACR: 5+/-0.7 mm\n", slice_thickness)
slice_position_accuracy = phantom.ACR_test_SlicePositionAccuracy(image.array(:,:,1));
sprintf("Slice Position Accuracy: %f mm vs ACR: -5<spa<5 mm\n", slice_position_accuracy)

[slice_uniformity, mean_pixel_intensity] = phantom.ACR_test_SliceUniformity(image.array(:,:,7));
sprintf("Percentage Integral Uniformity: %f %% vs ACR: >87.5 %% for 1.5T, >82 %% for 3T \n", slice_uniformity*100)
ghosting_ratio = phantom.ACR_test_PercentSignalGhosting(image.array(:,:,7));
sprintf("Ghosting Ratio: %f vs ACR: <0.025 ", ghosting_ratio)
if ghosting_ratio < 0.025; sprintf("PASS \n"); else sprintf("FAIL \n") endif

