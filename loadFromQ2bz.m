function img = loadFromQ2bz ( path )
    [pathstr, name, ext] = fileparts(path);
    if ~strcmp(ext,'.q2bz') && ~strcmp(ext,'.bz2')
      img = double(imread(path));
      img = img - min(img(:));
      return;
    end
    if isempty( pathstr )
      pathstr = '.';
    end
    tmpPath = [pathstr '/' name '.tmp'];   
    system(['bunzip2 -c ' path ' > ' tmpPath]);
    
    fid = fopen(tmpPath, 'r');
    fgetl(fid); % Skip magic number
    fgetl(fid); % Skip header
    
    % Read width and height
    arr = regexp(fgetl(fid),'\s+','split');
    width = str2num(arr{1});
    height = str2num(arr{2});
    
    % Skip max, but be careful not to read more than one new line after max.
    % The binary data could start with a value that is equivalent to a
    % new line.
    fscanf(fid,'%d',1);
    % This skips the character following max, i.e. the new line we are
    % expecting.
    fscanf(fid,'%c',1);
    
    % Read image to vector
    x = fread(fid,Inf,'double');
    
    img = zeros(width,height);
    for i=1:width
        for j=1:height
            img(j,i) = x(width*(j-1)+i);
        end;
    end;
    
    delete(tmpPath);
end