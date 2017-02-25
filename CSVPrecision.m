function [preclist] = CSVPrecision(path)
%% Load .csv file and .png file from under path
    % folder contain one image and the corresponding csv only
    % path = 'D:\2017\BumpFitting\SNR test\S21\50Frames_Results\B1';
    files = dir(fullfile(path,'*.csv'));
    images = dir(fullfile(path,'*.png'));
    files = {files.name};
    images = {images.name};

    for numframe = 1:size(images,2)
        filename = files{1,numframe};
        imagename = images{1,numframe};
        temp = csvread(fullfile(path,filename),1,0); %skip row #0
        img = imread(fullfile(path,imagename));

    %% filter the list to get only one set of sublattice
        % filter peaks with integer coordinates

        num=1;
        templist = zeros(1,2);
        for j = 1:size(temp,1)
            if mod(temp(j,1),1)<0.001 &&mod(temp(j,2),1)<0.0001
                templist(num) = j;
                num = num+1;
            end
        end
        temp(templist(:),:)=[];
        clear num;

        % filter peaks too close to border
        temp(temp(:,1)<5,:) = [];
        temp(temp(:,1)>145,:) = [];
        temp(temp(:,2)<5,:) = [];
        temp(temp(:,2)>145,:) = [];

        % filter to get one set of subset
        % click original point first, then one point along a direction and one
        % along b direction
        imshow(img,[]);hold on;
        scatter(temp(:,1)+1,temp(:,2)+1,30,'filled');
        [x_temp,y_temp] = ginput(3);
         close all;
        x_temp = double(x_temp);
        y_temp = double(y_temp);
        [x_temp(1),y_temp(1)] = findnearestpeak(x_temp(1),y_temp(1),temp);
        [x_temp(2),y_temp(2)] = findnearestpeak(x_temp(2),y_temp(2),temp);
        [x_temp(3),y_temp(3)] = findnearestpeak(x_temp(3),y_temp(3),temp);
        num=2;
        peaklist = [x_temp(1),y_temp(1);x_temp(2),y_temp(2)];
        separation_a = [x_temp(2)-x_temp(1);y_temp(2)-y_temp(1)];
        separation_b = [x_temp(3)-x_temp(1);y_temp(3)-y_temp(1)];

        % detection along a direction
        x_cor = x_temp(2);
        y_cor = y_temp(2);
        while ( (x_cor+separation_a(1))<145 && (y_cor+separation_a(2))<145 && num<6)
            x_cor = x_cor + separation_a(1);
            y_cor = y_cor + separation_a(2);
            [peaklist(num+1,1),peaklist(num+1,2)] = findnearestpeak(x_cor,y_cor,temp);
            num = num + 1;
        end

        % detection along b direction
        x_peaknum = size(peaklist,1);
        x_cor = x_temp(1);
        y_cor = y_temp(1);
        cycle = 0;
        while (x_cor+separation_b(1)<145 && y_cor+separation_b(2)<145 && cycle<5)
            for i = 1 : x_peaknum
                x_cor = peaklist(cycle * x_peaknum + i, 1);
                y_cor = peaklist(cycle * x_peaknum + i, 2);
                x_cor = x_cor + separation_b(1);
                y_cor = y_cor + separation_b(2);
                [peaklist(num+1,1),peaklist(num+1,2)] = findnearestpeak(x_cor, y_cor, temp);
                num = num + 1;
            end
            x_cor = peaklist(end,1);
            y_cor = peaklist(end,2);
            cycle = cycle + 1;
        end

%         imagesc(img); hold on;
%         scatter(peaklist(:,1)+1,peaklist(:,2)+1,30,'blue','filled');

        %% calculate precision along x and y direction with function Precision
        y_list = Precision(peaklist(:,1),peaklist(:,2),80,80,21,15,21.19,0); %calculate y separation
        x_list = Precision(peaklist(:,1),peaklist(:,2),10,-10,21,15,21.16,1); %calculate x separation
        x_prec = std(x_list(:,5));
        y_prec = std(y_list(:,5)); %precision in nm along x and y direction
        preclist(numframe,1)=x_prec;
        preclist(numframe,2)=y_prec;
    end

end

%% find peak position along list that is closest to (x_guess,y_guess), return the peak location in (x_loc,y_loc)
function [x_loc, y_loc] = findnearestpeak(x_guess, y_guess, list)
    for i = 1:size(list,1)
        list(i,3) = pdist([x_guess,y_guess;list(i,1),list(i,2)],'euclidean');
    end
    [~,Index] = min(list(:,3));
    x_loc = list(Index,1);
    y_loc = list(Index,2);
end