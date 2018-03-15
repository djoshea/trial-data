function imgout = convert2D(Conversion, img)
   % img and imgout is a x b x 3
   
    r = size(img, 1);
    c = size(img, 2);
    assert(size(img, 3) == 3);
    imgflat = reshape(img, [r*c 3]);
    imgoutflat = TrialDataUtilities.Color.convert(Conversion, imgflat);
    imgout = reshape(imgoutflat, [r c 3]);
end
