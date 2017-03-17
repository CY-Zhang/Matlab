% function used to find the distance between left-peak and right-peak for certain peak, opt_anlge=1 for horizontal and 0 for 
% vertical
function [separation_list] = SeparationAnalysis(X_location, Y_location, dist_max, px_size_x, px_size_y)

numpeaks = size(X_location,1);
separation_list = zeros(numpeaks,4); %1st column for x sep and 2nd for y sep

for i = 1:numpeaks
    
    list = zeros(numpeaks,4);
    list(:,1) = X_location;
    list(:,2) = Y_location;
    separation_list(i,3) = list(i,1);             %x coordinate
    separation_list(i,4) = list(i,2);             %y coordinate
    
    list(:,3) = abs(list(:,1) - list(i,1)); %x separation for all peaks
    list(:,4) = abs(list(:,2) - list(i,2));   %y separation for all peaks
    list(list(:,3)>dist_max,3) = 0; %filter x distance
    list(list(:,3)==0,:) = [];
    list(list(:,4)>dist_max,4) = 0; %filter y distance
    list(list(:,4)==0,:) = [];
    
    switch size(list,1)
        case 3 %corner
        separation_x = 0;
        separation_y = 0;
        case 5 %edge
            if sum(list(:,3))<sum(list(:,4)) %along y edge
                separation_x = 0;
                list_temp = sortrows(list,3);
                separation_y = pdist([list_temp(1,1:2);list_temp(2,1:2)],'euclidean');
                %separation_y = list_temp(1,4)+list_temp(2,4);
            else
                separation_y = 0;            %along x edge
                list_temp = sortrows(list,4);
                separation_x = pdist([list_temp(1,1:2);list_temp(2,1:2)],'euclidean');
                %separation_x = list_temp(1,3)+list_temp(2,3);
            end
        case 8 %center
            list_temp = sortrows(list,4);
            separation_x = list_temp(1,3)+list_temp(2,3);
            list_temp = sortrows(list,3);
            separation_y = list_temp(1,4)+list_temp(2,4);
    end
    
    
    separation_list(i,1) = separation_x.*px_size_x; %separation along x
    separation_list(i,2) = separation_y.*px_size_y; %separation along y

end

% plot separation along x direction
figure;
scatter(separation_list(separation_list(:,1)~=0,3),separation_list(separation_list(:,1)~=0,4),300,...
    separation_list(separation_list(:,1)~=0,1),'filled','s'); 
caxis([mean(separation_list(separation_list(:,1)~=0,1))-0.02 mean(separation_list(separation_list(:,1)~=0,1))+0.02]);
colorbar;

figure;
scatter(separation_list(separation_list(:,2)~=0,3),separation_list(separation_list(:,2)~=0,4),300,...
    separation_list(separation_list(:,2)~=0,2),'filled','s'); 
caxis([7.79 7.83]);
caxis([mean(separation_list(separation_list(:,2)~=0,2))-0.02 mean(separation_list(separation_list(:,2)~=0,2))+0.02]);
colorbar;

end