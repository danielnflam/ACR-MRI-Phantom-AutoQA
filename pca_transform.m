function img_transformed = pca_transform(X, img)
  # X is Nx2 array, with N observations.
  # img is the image to be transformed

  # PCA the points
  [coeff, score, ~] = pca(X);
  # score is the projection in the new coords --> use score = coeff * old coords
  coeff = [coeff'; 0 0];
  # Transformation of original image to img_transformed
  T = maketform("affine", coeff);
  img_transformed = imtransform(img,T);
endfunction

