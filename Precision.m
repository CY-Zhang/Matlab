function [distance_a]=Precision(X_location, Y_location, ang_max, ang_min, ...
    dist_max, dist_min, px_size, opt_angle)
%opt_anlge = 0, than use abs(angle)>ang_min as the only criteria for angle
%this option is helpful when angle is around 90 degrees, as it will take
%both positive and negative
%% Filter the X and Y location list, only peaks that is fitted would be left
% this part can be optimized
% num=1;
% for i = 1:size(X_location,1)
%     if mod(X_location(i),1)==0
%         templist(num)=i;
%         num = num+1;
%     end
% end
% X_location(templist(:))=[];
% Y_location(templist(:))=[];

%% Initialize variables
PI = 3.1415926;
numdumbbell = size(X_location,1);
num = 1;
dumbbell_list = zeros(numdumbbell,4);
dumbbell_list(:,1) = X_location;
dumbbell_list(:,2) = Y_location;

%%Loop through all peaks to find pairs with desired angle and separation
for i = 1:numdumbbell
    for j = 1:numdumbbell
        if(j~=i)
            x1 = X_location(i);
            y1 = Y_location(i);
            x2 = X_location(j);
            y2 = Y_location(j);
            if opt_angle == 1
                 if(x2>x1) %apply this criterial to avoid duplicate, atan the same if flip the order
                     %use x2>x1 when measuring separation along x direction
                    distance = sqrt((x2-x1)^2+(y2-y1)^2);
                    angle = atan((y2-y1)/(x2-x1))/PI*180;
                    if(opt_angle ==1)
                        if((angle>ang_min)&&(angle<ang_max)&&(distance>dist_min)&&(distance<dist_max))
                            dumbbell_list(i,3) = angle;
                            dumbbell_list(i,4) = distance;
                            num = num+1;
                        end
                    else
                        if((abs(angle)>ang_min)&&(distance>dist_min)&&(distance<dist_max))
                            dumbbell_list(i,3) = angle;
                            dumbbell_list(i,4) = distance;
                            num = num+1;
                        end
                    end
                 end
            end
            if opt_angle == 0
                if(y2>y1) 
                    distance = sqrt((x2-x1)^2+(y2-y1)^2);
                    angle = atan((y2-y1)/(x2-x1))/PI*180;
                    if(opt_angle ==1)
                        if((angle>ang_min)&&(angle<ang_max)&&(distance>dist_min)&&(distance<dist_max))
                            dumbbell_list(i,3) = angle;
                            dumbbell_list(i,4) = distance;
                            num = num+1;
                        end
                    else
                        if((abs(angle)>ang_min)&&(distance>dist_min)&&(distance<dist_max))
                            dumbbell_list(i,3) = angle;
                            dumbbell_list(i,4) = distance;
                            num = num+1;
                        end
                    end
                 end
            end
        end
    end
end
list_refined = dumbbell_list((dumbbell_list(:,3)~=0),:);
distance_a(:,1) = list_refined(:,1);%x coordinates
distance_a(:,2) = list_refined(:,2);%y coordinates
distance_a(:,3) = list_refined(:,3);%angle list
distance_a(:,4) = list_refined(:,4);%distance in px
distance_a(:,5) = distance_a(:,4).*px_size; %distance in nm

Figure1 = figure('position',[100,100,450,300]);
scatter(distance_a(:,1)+1,distance_a(:,2)+1,300,distance_a(:,5),'filled','s');
%caxis(ca_nm);
colorbar;
%ax = gca;
set(gca,'Ydir','reverse','FontSize',18);
axis equal;

end