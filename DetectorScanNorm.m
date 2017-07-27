%%Do Normalization based on Detector scan image
function HAADF_norm = DetectorScanNorm(BeamOn, BeamOff, HAADF, nbins)
%%First average over the scan map to get BeamLevel and DarkLevel
addpath 'D:\MatlabCode'
load DetectorResponse.mat;

%% calculate beam level with threshold
Threshold = 5000;
mask = heaviside(BeamOn-Threshold);
Detector_mask = BeamOn.*mask;
Beam_level = mean(Detector_mask(Detector_mask~=0));

%% calculate dark level from center disk with 20px radius
cx=1040; cy=1021;
ix=2048; iy=2048;
r=20;
[x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
background_mask=((x.^2+y.^2)<=r^2);
background_mask = background_mask .* BeamOn;
imagesc(background_mask);
Dark_level = mean(background_mask(background_mask~=0));

%% normalize HAADF image and print results
HAADF_norm = (HAADF-Dark_level)/Beam_level;
sensmap = BeamOn./Beam_level;
subplot(1,2,1);imagesc(HAADF_norm);axis equal off; colorbar;
subplot(1,2,2);imagesc(sensmap);axis equal off; colorbar;
end