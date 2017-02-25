function img = loadFromQ3 ( path )
    [pathstr, name, ext] = fileparts(path);
    if isempty( pathstr )
      pathstr = '.';
    end
    fid = fopen(path, 'r');
    fgetl(fid); % Skip magic number
    fgetl(fid); % Skip header
    
    % Read width and height
    arr = regexp(fgetl(fid),'\s+','split');
    width = str2num(arr{1})
    height = str2num(arr{2})
    depth = str2num(arr{3})
    
    % Skip max, but be careful not to read more than one new line after max.
    % The binary data could start with a value that is equivalent to a
    % new line.
    fscanf(fid,'%d',1);
    % This skips the character following max, i.e. the new line we are
    % expecting.
    fscanf(fid,'%c',1);
    
    % Read image to vector
    x = fread(fid,Inf,'double');
    img = reshape ( x, width, height, depth );
end