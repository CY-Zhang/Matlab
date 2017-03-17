function [rect_list, center_list] = CenterAnalysis(Sr_X, Sr_Y, Ti_X, Ti_Y, HAADF, mode)
x_sep = 19;
y_sep = 19;
delta = 3;
col_y = 9;
col_x = 11;
px_sizex = 0.2101;
px_sizey = 0.2131;

if mode == 2
    num_atom = size(Ti_X,1);
    rect_list = zeros(num_atom,10);
    center_list = zeros(num_atom,7);
    for i = 1:num_atom
        % initialize coordinates with initial guess
        x1 = Ti_X(i);
        y1 = Ti_Y(i);
        x2 = x1 + x_sep;
        y2 = y1;
        x3 = x1 + x_sep;
        y3 = y1 + y_sep;
        x4 = x1;
        y4 = y1 + y_sep;
        % break if atom1 is on top or right border
        if (x3 > max(Ti_X)+5 || y3 > max(Ti_Y)+5) %no atom on top or right
            continue
        end
        % refine coordinates of atom2 - atom4
        x2 = Ti_X(logical((abs(Ti_X(:)-x2)<delta) .* (abs(Ti_Y(:)-y2)<delta)));
        y2 = Ti_Y(logical((abs(Ti_X(:)-x2)<delta) .* (abs(Ti_Y(:)-y2)<delta)));
        x3 = Ti_X(logical((abs(Ti_X(:)-x3)<delta) .* (abs(Ti_Y(:)-y3)<delta)));
        y3 = Ti_Y(logical((abs(Ti_X(:)-x3)<delta) .* (abs(Ti_Y(:)-y3)<delta)));
        x4 = Ti_X(logical((abs(Ti_X(:)-x4)<delta) .* (abs(Ti_Y(:)-y4)<delta)));
        y4 = Ti_Y(logical((abs(Ti_X(:)-x4)<delta) .* (abs(Ti_Y(:)-y4)<delta)));

        rect_list(i,1) = x1;
        rect_list(i,2) = x2;
        rect_list(i,3) = x3;
        rect_list(i,4) = x4;
        rect_list(i,5) = y1;
        rect_list(i,6) = y2;
        rect_list(i,7) = y3;
        rect_list(i,8) = y4;
        rect_list(i,9) = (x3+x4+x1+x2)/4;
        rect_list(i,10) = (y3+y4+y1+y2)/4;
        
        if sum(logical((Sr_X < x3) .* (Sr_X > x1).* (Sr_Y <y3) .* (Sr_Y > y1)))
            center_list(i,1) = Sr_X(logical((Sr_X < x3) .* (Sr_X > x1).* (Sr_Y <y3) .* (Sr_Y > y1)));
            center_list(i,2) = Sr_Y(logical((Sr_X < x3) .* (Sr_X > x1).* (Sr_Y <y3) .* (Sr_Y > y1)));
        else
            continue
        end
        center_list(i,3) = pdist([x1*px_sizex,y1*px_sizey;center_list(i,1)*px_sizex,center_list(i,2)*px_sizey],'euclidean');
        center_list(i,4) = pdist([x2*px_sizex,y2*px_sizey;center_list(i,1)*px_sizex,center_list(i,2)*px_sizey],'euclidean');
        center_list(i,5) = pdist([x3*px_sizex,y3*px_sizey;center_list(i,1)*px_sizex,center_list(i,2)*px_sizey],'euclidean');
        center_list(i,6) = pdist([x4*px_sizex,y4*px_sizey;center_list(i,1)*px_sizex,center_list(i,2)*px_sizey],'euclidean');
        center_list(i,7) = (center_list(i,3)+center_list(i,4)+center_list(i,5)+center_list(i,6))/4;
        % how to define the center if the rectangular is not perfectly
        % vertical?
    end
    
    
    center_list = center_list(center_list(:,1)~=0,:);
    
    figure;
    scatter(center_list(:,1),center_list(:,2),300,center_list(:,7),'filled');
    colorbar;
    
    [x_temp,index] = sort(center_list(:,1));
    y_temp = center_list(index,2);
    z_temp = center_list(index,7);

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
    colormap hot;
    axis equal off;
end



if mode == 1 %calculate Ti offcenter and use quiver to plot the offcenter amount
    num_atom = size(Sr_X,1);
    rect_list = zeros(num_atom,10);
    center_list = zeros(num_atom,2);
    for i = 1:num_atom
        % initialize coordinates with initial guess
        x1 = Sr_X(i);
        y1 = Sr_Y(i);
        x2 = x1 + x_sep;
        y2 = y1;
        x3 = x1 + x_sep;
        y3 = y1 + y_sep;
        x4 = x1;
        y4 = y1 + y_sep;
        % break if atom1 is on top or right border
        if (x3 > max(Sr_X)+5 || y3 > max(Sr_Y)+5) %no atom on top or right
            continue
        end
        % refine coordinates of atom2 - atom4
        x2 = Sr_X(logical((abs(Sr_X(:)-x2)<delta) .* (abs(Sr_Y(:)-y2)<delta)));
        y2 = Sr_Y(logical((abs(Sr_X(:)-x2)<delta) .* (abs(Sr_Y(:)-y2)<delta)));
        x3 = Sr_X(logical((abs(Sr_X(:)-x3)<delta) .* (abs(Sr_Y(:)-y3)<delta)));
        y3 = Sr_Y(logical((abs(Sr_X(:)-x3)<delta) .* (abs(Sr_Y(:)-y3)<delta)));
        x4 = Sr_X(logical((abs(Sr_X(:)-x4)<delta) .* (abs(Sr_Y(:)-y4)<delta)));
        y4 = Sr_Y(logical((abs(Sr_X(:)-x4)<delta) .* (abs(Sr_Y(:)-y4)<delta)));

        rect_list(i,1) = x1;
        rect_list(i,2) = x2;
        rect_list(i,3) = x3;
        rect_list(i,4) = x4;
        rect_list(i,5) = y1;
        rect_list(i,6) = y2;
        rect_list(i,7) = y3;
        rect_list(i,8) = y4;
        rect_list(i,9) = (x3+x4+x1+x2)/4;
        rect_list(i,10) = (y3+y4+y1+y2)/4;

        center_list(i,1) = Ti_X(logical((Ti_X < x3) .* (Ti_X > x1).* (Ti_Y <y3) .* (Ti_Y > y1)));
        center_list(i,2) = Ti_Y(logical((Ti_X < x3) .* (Ti_X > x1).* (Ti_Y <y3) .* (Ti_Y > y1)));
        % how to define the center if the rectangular is not perfectly
        % vertical?
    end

    center_list(center_list(:,1)==0,:)=[];
    rect_list(rect_list(:,1)==0,:)=[];
    displacement = center_list - rect_list(:,9:10);
    x = reshape(center_list(:,1),[num_col,num_row]);
    x = fliplr(x');
    y = reshape(center_list(:,2),[num_col,num_row]);
    y = fliplr(y');
    u = reshape(displacement(:,1),[num_col,num_row]);
    u = fliplr(u');
    v = reshape(displacement(:,2),[num_col,num_row]);
    v = fliplr(v');
    % quiver(x,y,u,v);

    figure;
    imshow(HAADF',[]);
    hold on;
    ax1 = quiver(x+1,y+1,u,v);
    set(ax1,'LineWidth',2);
    set(gca,'YDir','reverse');
    axis xy;
end
end