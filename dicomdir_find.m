function OUT = dicomdir_find(DICOMDIR)
  # find files in DICOM folder
  count = 0;
  OUT = {};
  for i = 3:length(DICOMDIR)
    if ~DICOMDIR(i).isdir && DICOMDIR(i).name != "VERSION"
      count = count+1;
      OUT{count} = DICOMDIR(i).name;
    endif
  end

endfunction
