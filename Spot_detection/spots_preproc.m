function [dataCells, dataSpots, timepoints, xpoints, xpos_idx, direction, ypoints, y_heights, y_groups, fovNum,...
    nCellsFOV, nCells, spotFrac, spotFracFOV, zeroSpot, oneSpot, twoPlusSpot, sosFOV, cellLenFOV, cmap] =...
    spots_preproc(Bacmman_folder, dataset_name, prefix, files_folder, bf_keyword)

% Import data
dataCells = readtable([Bacmman_folder dataset_name '/' dataset_name '_0.csv'], 'TreatAsEmpty', 'NA');
dataSpots = readtable([Bacmman_folder dataset_name '/' dataset_name '_1.csv'], 'TreatAsEmpty', 'NA');

% Sort datasets chronologically
for i = 1:height(dataCells)
    name = dataCells.Position(i);
    if ~isempty(prefix)
        name = char(name);
        idx = name(strfind(name, prefix)+length(prefix):end);
        dataCells.TrueIdx(i) = str2double(idx);
    else
        dataCells.TrueIdx(i) = name;
    end
end
dataCells = sortrows(dataCells, 'TrueIdx', 'ascend');
for i = 1:height(dataSpots)
    name = dataSpots.Position(i);
    if ~isempty(prefix)
        name = char(name);
        idx = name(strfind(name, prefix)+length(prefix):end);
        dataSpots.TrueIdx(i) = str2double(idx);
    else
        dataSpots.TrueIdx(i) = name;
    end
end
dataSpots = sortrows(dataSpots, 'TrueIdx', 'ascend');

% Fetch timestamps and stage position for each FOV
cd(files_folder)
listing = dir(['*' bf_keyword '*']);
uCellIdx = unique(dataCells.TrueIdx, 'sorted');
for i = 1:length(uCellIdx)
    bfInfo = imfinfo(listing(i).name);
    dataCells.Timestamp(dataCells.TrueIdx == uCellIdx(i)) = datetime(bfInfo(1).DateTime, 'InputFormat', 'yyyyMMdd HH:mm:ss.SSS');
    ASI_pos = strfind(bfInfo(1).ImageDescription,'ASI X');
    substr = bfInfo(1).ImageDescription(ASI_pos:ASI_pos+100);
    val_pos = strfind(substr,'value="');
    end_pos = strfind(substr,'"/>');
    asi_x = str2double(substr(val_pos(1)+7:end_pos(1)-1));
    asi_y = str2double(substr(val_pos(2)+7:end_pos(2)-1));
    dataCells.asi_x(dataCells.TrueIdx == uCellIdx(i)) = asi_x;
    dataCells.asi_y(dataCells.TrueIdx == uCellIdx(i)) = asi_y;
end
dataCells.Time = dataCells.Timestamp - dataCells.Timestamp(1);
timepoints = unique(dataCells.Time);
xpoints = grpstats(dataCells.asi_x, dataCells.TrueIdx, 'mean');
xpoints_original = xpoints;
direction = [0; diff(xpoints_original) > 0]; % 0: left-to-right; 1: right-to-left
[xpoints, xpos_idx] = sort(xpoints);
direction = direction(xpos_idx);
ypoints = round(grpstats(dataCells.asi_y, dataCells.TrueIdx, 'mean'));
y_groups = findgroups(round(ypoints));
ypoints = ypoints(xpos_idx);
y_heights = unique(round(ypoints));


% Measures
fovNum = unique(dataCells.TrueIdx);
nCellsFOV = grpstats(dataCells.Idx, dataCells.TrueIdx, 'max') +1;
nCells = height(dataCells);
spotFrac = sum(dataCells.SpotCount > 0)/nCells;
spotFracFOV = grpstats(dataCells.SpotCount > 0, dataCells.TrueIdx, 'sum')./nCellsFOV;
zeroSpot = grpstats(dataCells.SpotCount == 0, dataCells.TrueIdx, 'sum')./nCellsFOV;
oneSpot = grpstats(dataCells.SpotCount == 1, dataCells.TrueIdx, 'sum')./nCellsFOV;
twoPlusSpot = grpstats(dataCells.SpotCount >= 2, dataCells.TrueIdx, 'sum')./nCellsFOV;

dataCells.normSOS = dataCells.MeanSOS./dataCells.SpineLength;
if ismember('MeanSOS', dataCells.Properties.VariableNames)
    sosFOV = grpstats(dataCells.normSOS, dataCells.TrueIdx, 'mean');
end
cellLenFOV = grpstats(dataCells.SpineLength, dataCells.TrueIdx, 'mean');


% Normalise spot positions in the cell
% SpineCurvilinearCoord: [0:cellLength]
% SpineRadialCoord: distance from spine (pixels)

dataSpots.normXpos = dataSpots.SpineCurvilinearCoord./dataSpots.SpineLength -0.5; % Normalise to cell length & centre
dataSpots.normYpos = dataSpots.SpineRadialCoord./dataSpots.SpineRadius; % Normalise to cell radius

cmap = rand(length(y_heights), 3);

end