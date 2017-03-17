% find sourcesize and sample thickness that would maximize the xcorrelation
% with given exp and simulation stack
function [SS, thickness] = StackSSCalibrate(exp)
 % parameters
 plotresult = 0;
 min_uc = 51;   % min and max allowed thickness
 max_uc = 51;
 plotconvergence = 1;
 
% load in the simulation stack
 addpath('D:\2017\STO_SRO');
 sim_stack = IBWread('HAADF_stack.ibw');
 sim_stack = sim_stack.y;
 xcorrlist = zeros(size(sim_stack,3),1);
 SSlist = xcorrlist .* 0;

 
 for layer = min_uc : max_uc
     [SSlist(layer), xcorrlist(layer)] = SourceSizeCalibrate(exp,sim_stack(:,:,layer),plotconvergence);
 end
 
 [~, thickness] = min(xcorrlist);
 SS = SSlist(thickness);
 
 if plotresult
    figure; 
    subplot(2,1,1);
    plot(SSlist,'LineWidth',2,'Color',[0    0.4470    0.7410]);
    xlabel('Number of Multislice Layers');
    ylabel('Optimized Source Size (Angstrom)');
%     ylim([0 1]);
    hold on;
    plot([thickness thickness],[0 Inf],'LineWidth',0.5,'Color','b');
    
    subplot(2,1,2);
    plot(-xcorrlist,'LineWidth',2,'Color',[0.8500    0.3250    0.0980]);
    xlabel('Number of Multislice Layers');
    ylabel('Averaged XCorrelation for All Peaks');
%     ylim([0 1]);
    hold on;
    plot([thickness thickness],[0 Inf],'LineWidth',0.5,'Color','b');

 end

 thickness = thickness/2; % Multislice layers -> uc (for STO)
 
 fprintf('The optimized source size is %.3f Angstrom\n',SS);
 fprintf('Sample thickness is %d layers of unit cell\n',thickness);
end

% find sourcesize that would maximize the xcorrelation with given exp and
% simulation image
function [SourceSize, avg_xcorr] = SourceSizeCalibrate(exp, sim, plotoption)

exp_pxsize = 21.16;
sim_pxsize = 9.7625;

ratio = sim_pxsize/exp_pxsize; %sim_pxsize/exp_pxsize
sim_resize = imresize(sim,ratio);
C = normxcorr2(sim_resize,exp);
threshold = heaviside(C-0.75);

peak = imregionalmax(C.*threshold);
peaklist = zeros(sum(peak(:)),2);
[peaklist(:,1), peaklist(:,2)] = find(peak==1);
peaklist(:,1) = peaklist(:,1) - size(sim_resize, 1)+1;
peaklist(:,2) = peaklist(:,2) - size(sim_resize, 2)+1;


ss_size = 1.6; % initial guess of standard deviation in px

if plotoption
    options = optimset('PlotFcns',@optimplotfval,'display','iter'); 
else
    options = optimset('display','off');
end

func = @(ss) averagexcorr(ss,exp,sim_resize,peaklist(:,1),peaklist(:,2));
[SourceSize, avg_xcorr] = fminsearch(func, ss_size, options); % source size in pixel
SourceSize = SourceSize*2.35482*21.16/100; %relationship between FWHM and std, then into angstrom
% SourceSize = SourceSize*2.35482*9.7625/100;

% fprintf('The optimized source size is %.3f Angstrom\n',SourceSize);


end

function avg_corr = averagexcorr(ss,exp,sim,row,col)
    %sim_ss = imgaussfilt(sim, [ss ss]); %replace imgaussfilt with more stable way
    if ss < 0
        avg_corr = Inf;
    else
%         sim_fft = fft2(sim);
%         sim_fft = fftshift(sim_fft);
%         sim_fft = sim_fft .* fspecial('gaussian',size(sim_fft),ss);
%         sim_ss = abs(ifft2(sim_fft));
        sim_ss = imgaussfilt(sim, ss);
        max_list = row.*0;
        for i = 1:size(row,1)
             temp = exp(row(i):row(i)+size(sim,1)-1,col(i):col(i)+size(sim,2)-1);
%             max_list(i) = -ssim(temp, sim_ss);
%             diff = (sim_ss - temp).^2;
%             max_list(i) = sum(diff(:));
            C = xcorr2(sim_ss, temp);
            max_list(i) = -max(C(:));
        end
        avg_corr = mean(max_list);
    end
end

