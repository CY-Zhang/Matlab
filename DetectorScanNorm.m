%%Do Normalization based on Detector scan image
function HAADF_norm = DetectorScanNorm(BeamOn, BeamOff, HAADF, nbins)
%%First average over the scan map to get BeamLevel and DarkLevel
addpath 'D:\MatlabCode'
load DetectorResponse.mat;
Detector_hist = histogram(BeamOn,nbins);
histc = histcounts(BeamOn,nbins);
histc_diff = diff(histc);
[~, MinIdx] = min(abs(histc_diff));
Threshold = (MinIdx+1)*Detector_hist.BinWidth;
Threshold = 10000;
mask = heaviside(BeamOn-Threshold);
Detector_mask = BeamOn.*mask;
Background_mask = BeamOn.*(1-mask); %BG calculated in this way is
%significantly higher than expected
Beam_level = mean(BeamOn(mask~=0));
Dark_level = mean(BeamOn(mask==0));
Dark_level = mean(mean(BeamOff));
HAADF_norm = (HAADF-Dark_level)/Beam_level;
sensmap = BeamOn./Beam_level;
subplot(1,2,1);imagesc(HAADF_norm);axis equal off; colorbar;
subplot(1,2,2);imagesc(sensmap);axis equal off; colorbar;
end