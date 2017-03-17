function Int = PeakAnalysis

%% Atomic Counting
clear all
plotresult = 1;
col_x = 11;
col_y = 9;
path = 'G:\030417\S1\HAADF';

% addpath('C:\Users\Ryan\Dropbox\STO_SRO\S3\HAADF');
addpath('G:\030417\StandardHAADF\Sr');
% addpath('C:\Users\Chenyu\Dropbox\STO_SRO\S3\HAADF');%simulation path
% addpath('D:\2017\STO_SRO\030417\S7');% exp path

% load simulation fitting parameters, should be the same for all STO series
A_sim = IBWread('A_stack.ibw');
A_sim = A_sim.y;

Xw_sim = IBWread('xW_stack.ibw');
Xw_sim = Xw_sim.y;
Xw_sim = Xw_sim.*9.76; %convert to physical unit pm

Yw_sim = IBWread('yW_stack.ibw');
Yw_sim = Yw_sim.y;
Yw_sim = Yw_sim.*9.76;

Cor_sim = IBWread('cor_stack.ibw');
Cor_sim = Cor_sim.y;

z0_sim = IBWread('z0_stack.ibw');
z0_sim = z0_sim.y;

% load experiment fitting parameters
w1 = fullfile(path,'A.ibw');
w1 = IBWread(w1);
w1 = w1.y;

w2 = fullfile(path,'x0.ibw');
w2 = IBWread(w2);
w2 = w2.y;

w3 = fullfile(path,'xW.ibw');
w3 = IBWread(w3);
w3 = w3.y;
w3_pm = w3.*21.16; %convert to physical unit pm

w4 = fullfile(path,'y0.ibw');
w4 = IBWread(w4);
w4 = w4.y;

w5 = fullfile(path,'yW.ibw');
w5 = IBWread(w5);
w5 = w5.y;
w5_pm = w5.*21.19;

w6 = fullfile(path,'cor.ibw');
w6 = IBWread(w6);
w6 = w6.y;

w0 = fullfile(path,'z0.ibw');
w0 = IBWread(w0);
w0 = w0.y;

fclose('all');
% windowsize = 50; %window size in pm for Sr
% windowsize = 50; %window size in pm

% intensity following LeBeau's scheme, works well with half window size 112pm
% Int_simulation = A_sim.*0;
% for i = 1:size(Int_simulation,2)
%     Int_simulation(i) = GaussIntegrate(windowsize,z0_sim(i), A_sim(i),Xw_sim(i),Yw_sim(i),Cor_sim(i));
% end
% Int_simulation = Int_simulation';
% 
% 
% Int = zeros(size(w1,1),4);
% for i = 1:size(w1,1)
%     Int(i,3) = GaussIntegrate(windowsize, w0(i), w1(i), w3_pm(i), w5_pm(i), w6(i));
%     temp = abs(Int_simulation - Int(i,3));
%     [~,index] = min(temp);
%     Int(i,4) = index; % number of multislice layer
%     Int(i,4) = index/2;
%     %Int(i,4) = floor(Int(i,4)/2);
% end
% clear temp index

% intensity following Jie's scheme
windowsize = 40;
% windowsize = 70;
Int_simulation = A_sim.*0;
for i = 1:size(Int_simulation,2)
    Int_simulation(i) = z0_sim(i)*windowsize*windowsize + 2*A_sim(i)*pi*Xw_sim(i)*Yw_sim(i)*sqrt(1-Cor_sim(i)^2);
end

Int = zeros(size(w1,1),3);
for i = 1:size(w1,1)
    Int(i,3) = w0(i)*windowsize*windowsize + 2*w1(i)*pi*w3_pm(i)*w5_pm(i)*sqrt(1-w6(i)^2);
    temp = abs(Int_simulation - Int(i,3));
    [~,index] = min(temp);
    Int(i,4) = index; % number of multislice layer
    Int(i,4) = index/2;
    %Int(i,4) = floor(Int(i,4)/2);
end
clear temp index

Int(:,1) = w2;
Int(:,2) = w4;
[Int(:,1),index] = sort(Int(:,1));
Int(:,2) = Int(index,2);
Int(:,3) = Int(index,3);
Int(:,4) = Int(index,4); %sort whole sheet by first column, i.e. x cord

if plotresult
    % simple plot way using scatter
    figure;
    X = Int(:,1);
    Y = Int(:,2);
    Z = Int(:,4);
    scatter(X,Y,300,Z,'s','filled');
    colorbar;
    % plot result with 3D bar plot, need to first sort all x and y data
    figure;
    X = reshape(Int(:,1),[col_y col_x]);
    Y = reshape(Int(:,2),[col_y col_x]);
    Z = reshape(Int(:,4),[col_y col_x]); 
    for j = 1:col_x
        temp = Y(:,j);
        [temp,index] = sortrows(temp);
        Y(:,j) = temp;
        temp = X(:,j);
        temp = temp(index);
        X(:,j) = temp;
        temp = Z(:,j);
        temp = temp(index);
        Z(:,j) = temp;
    end
    %Z = flipud(Z);  % now Z(1,1) corresponds to the bottom left corner Sr
    [r,c] = size(Z);                    %# Size of Z
    Y = 1:col_y;                        %# The positions of bars along the y axis
    C = mat2cell(kron(Z,ones(6,4)),6*r,4.*ones(1,c)).';  %'# Color data for Z
    
    hBar = bar3(Y,Z);           %# Create the bar graph
    set(gca,'Xdir','reverse');
    set(hBar,{'CData'},C);      %# Add the color data
    view(0,90);  %# Change the camera view: (0,90) for a 2D projection view
    view(157,51);
    colorbar;                   %# Add the color bar
end

% %% separation analysis along y direction
% temp = Precision(w2,w4,80,80,21,16,0.2119,0);
% for i = 1:size(temp,1)
%     coor = [temp(i,1) temp(i,2)];
%     intensity = Int(logical((Int(:,1)==temp(i,1)).* (Int(:,2)==temp(i,2))),3);
%     % what should be done here is if angle > 0, use -cos; if angle < 0, use
%     % cos(pi-angle)
%     coor_neighbor = [-temp(i,4)*cos(temp(i,3)*pi/180) abs(temp(i,4)*sin(temp(i,3)*pi/180))];
%     coor_neighbor = coor_neighbor + coor;
%     int_neighbor = Int(logical((abs(Int(:,1)-coor_neighbor(1))<1) .* (abs(Int(:,2)-coor_neighbor(2))<1)),3);
%     int_change = intensity - int_neighbor;
%     temp(i,6) = int_change;
% end
% 
% figure;
% scatter(temp(:,5),temp(:,6),'filled');
% xlabel('Separation along Y (Angstrom)');
% ylabel('Peak Intensity Change');
% 
% %% separation analysis along x direction
% temp = Precision(w2,w4,10,-10,21,16,0.2116,1);
% for i = 1:size(temp,1)
%     coor = [temp(i,1) temp(i,2)];
%     intensity = Int(logical((Int(:,1)==temp(i,1)).* (Int(:,2)==temp(i,2))),3);
%     coor_neighbor = [temp(i,4)*cos(temp(i,3)*pi/180) temp(i,4)*sin(temp(i,3)*pi/180)];
%     coor_neighbor = coor_neighbor + coor;
%     int_neighbor = Int(logical((abs(Int(:,1)-coor_neighbor(1))<1) .* (abs(Int(:,2)-coor_neighbor(2))<1)),3);
%     int_change = intensity - int_neighbor;
%     temp(i,6) = int_change;
% end
% 
% figure;
% scatter(temp(:,5),temp(:,6),'filled');
% xlabel('Separation along X (Angstrom)');
% ylabel('Peak Intensity Change');
end