function [OUT, OUT_meta] = dicomdir_read(DICOMDIR_PATH)
  # Reads out a DICOM directory containing multiple DICOM files

  # find files in DICOM folder
  count = 0;
  FILES = {};
  DICOMDIR = dir(DICOMDIR_PATH);

  for i = 3:length(DICOMDIR)
    if ~DICOMDIR(i).isdir && ~strcmp(DICOMDIR(i).name, "VERSION")
      count = count+1;
      FILES{count} = DICOMDIR(i).name;
    endif
  endfor

  # Order the slices in the correct order based on DICOM metadata SliceLocation
  sliceLoc = zeros(length(FILES),1);
  for i = 1:length(FILES)
    dicomfile_path = horzcat(DICOMDIR_PATH,filesep(),FILES{i});
    dicom_metadata = dicominfo(dicomfile_path);
    sliceLoc(i) = dicom_metadata.SliceLocation;
  endfor
  [sortedSliceLoc,idx] = sort(sliceLoc, "ascend");

  # Read out data into a 3D matrix
  OUT = zeros(dicom_metadata.Rows, dicom_metadata.Columns, length(FILES));
  count = 0;
  for i = idx'
    count = count+1;
    dicomfile_path = horzcat(DICOMDIR_PATH,filesep(),FILES{i});
    OUT(:,:,count) = dicomread(dicomfile_path);
  endfor

  OUT_meta.dcmMetadata = dicom_metadata;
  OUT_meta.sliceLocations = sortedSliceLoc;

  endfunction
