function varargout = get_features(folder, filename, position, varargin)

if length(varargin) < 1
    error('get_features requires at least a path, filename, image number and 1 feature name to be extracted')
end

if nargout ~= length(varargin)
    error('Number of output variables for get_features must be equal to the number of extracted features')
end

info = h5info([folder filename]);
dataset = info.Groups(1).Groups.Groups(position).Name;

for i = 1:length(varargin)
    
    feature = varargin{i};
    
    try
        data = h5read([folder filename],...
        [dataset '/' feature]);
    catch
        error(['Could not find feature "' feature '" for position ' num2str(position)])
    end
    
    outData = struct;
    outData.Connectivity = 8;
    [sz1, sz2, ~] = size(data); % For compatibility: "size(data, [1, 2])" introduced in R2019b
    outData.ImageSize = [sz1 sz2];
    outData.NumObjects = max(data, [], 'all');
    outData.PixelIdxList = cell(1, outData.NumObjects);
    
    for j=1:outData.NumObjects
        outData.PixelIdxList{j} = find(data == j);
    end
    
    varargout{i} = outData;

    
end


end