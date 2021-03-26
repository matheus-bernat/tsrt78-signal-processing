function x=resample(x,U,D)
%Resampling a factor U/D
   x=upsample(x,U);
   x=downsample(x,D);

