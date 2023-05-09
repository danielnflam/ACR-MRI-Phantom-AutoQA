classdef Image3D
  properties
    array # image array
    pixelSpacing # mm
    imagePlaneOrigin # mm
  endproperties

  methods
    function image = Image3D(d3d, sliceLocations, PixelSpacing, ImagePlaneOrigin)
      image.array = d3d;
      image.pixelSpacing = [PixelSpacing(1), PixelSpacing(2), abs(sliceLocations(2)-sliceLocations(1))];
      image.imagePlaneOrigin = ImagePlaneOrigin;
    endfunction

    function display(image, slice_num)
      imagesc(image.array(:,:,slice_num))
    endfunction

  endmethods
endclassdef
