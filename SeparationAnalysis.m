% function used to find the distance between left-peak and right-peak for certain peak, opt_anlge=1 for horizontal and 0 for 
% vertical
function SeparationAnalysis(X_location, Y_location, dist_max, px_size)

numpeaks = size(X_location,1);
separation_list = zeros(numpeaks,2); %1st column for x sep and 2nd for y sep

for i = 1:numpeaks
    
    list = zeros(numpeaks,4);
    list(:,1) = X_location;
    list(:,2) = Y_location;
    list(:,3) = abs(list(:,1) - list(i,1)); %x separation for all peaks
    list(:,4) = abs(list(:,2) - list(i,2)); %y separation for all peaks
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
                separation_y = list_temp(1,4)+list_temp(2,4);
            else
                separation_y = 0;            %along x edge
                list_temp = sortrows(list,4);
                separation_x = list_temp(1,3)+list_temp(2,3);
            end
        case 8
            list_temp = sortrows(list,4);
            separation_x = list_temp(1,3)+list_temp(2,3);
            list_temp = sortrows(list,3);
            separation_y = list_temp(1,4)+list_temp(2,4);
    end
    
    
    separation_list(i,1) = separation_x.*px_size;
    separation_list(i,2) = separation_y.*px_size;
end

end