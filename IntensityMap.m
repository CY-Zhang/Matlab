function [list] = IntensityMap(HAADF,x,y)
numatoms = size(x,1);
px_sizex = 21.01;
px_sizey = 21.31;
int_range_x = 390.5/4/px_sizex;%inteagration range along x and y
int_range_y = 390.5/4/px_sizey;
radius = int_range_x;
col_y = 8;
col_x = 5;

list = zeros(numatoms,3);
list(:,1) = x;
list(:,2) = y;
for i = 1:numatoms
    centre=[x(i) y(i)]; %center = [row col]
    Disk = fspecial('disk',radius)==0;
    mask = zeros(size(HAADF));
    if (centre(1)<radius || centre(1)-radius+size(Disk,1)-1 > size(mask,1) || ...
            centre(2) < radius || centre(2)-radius+size(Disk,2)-1 > size(mask,2))
        continue
    end
    mask(centre(1)-radius:centre(1)-radius+size(Disk,1)-1, centre(2)-radius:centre(2)-radius+size(Disk,2)-1)=double(~Disk);
    mask = mask(1:size(HAADF,1),1:size(HAADF,2));
    HAADF_masked = HAADF.*mask;
    list(i,3) = sum(HAADF_masked(:))/sum(mask(:));
end
list(list(:,3)==0,:)=[];

figure;
scatter(list(:,1),list(:,2),300,...
    list(:,3),'filled','s'); 
colorbar;
colormap hot;

[x_temp,index] = sort(list(:,1));
y_temp = list(index,2);
z_temp = list(index,3);

x_temp = reshape(x_temp,[col_y col_x]);
y_temp = reshape(y_temp,[col_y col_x]);
z_temp = reshape(z_temp,[col_y col_x]); 
for j = 1:col_x
    temp = y_temp(:,j);
    [temp,index] = sortrows(temp);
    y_temp(:,j) = temp;
    temp = x_temp(:,j);
    temp = temp(index);
    x_temp(:,j) = temp;
    temp = z_temp(:,j);
    temp = temp(index);
    z_temp(:,j) = temp;
end
% temp = reshape(distance_a(:,5),[7 10]);
figure;
imagesc(flipud(z_temp));
colorbar;
colormap gray;
axis equal off;
end