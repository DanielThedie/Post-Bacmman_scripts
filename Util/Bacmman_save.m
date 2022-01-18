% Bacmman_save
% Use this script to save Bacmman pre-processing or segmentation data to another local folder
% Data can be re-loaded to the Bacmman folder when necessary
% 
% CAUTION: loading will erase any existing Bacmman pre-processing/segmentation!

bacmman_folder = 'D:\Daniel\BACMMAN\Timelapse\210225_no_gRNA';
save_folder = 'D:\Daniel\Bacmman_save';

save_preprocessing = 1;
save_segmentation = 1;

load_preprocessing = 0;
load_segmentation = 0;

%% End of input

% List folders to be processed
ls = dir([bacmman_folder '/Output']);

% Safety checks
checks_passed = 'y';
if (save_preprocessing && load_preprocessing) || (save_segmentation && load_segmentation)
    error('Cannot save and load at the same time - check your input')    
end
if (load_preprocessing || load_segmentation) && length(dir([save_folder '/Output'])) < 3
    error(['No saved data found in ' save_folder]);
elseif (load_preprocessing || load_segmentation) && length(dir([save_folder '/Output']))-2 < length(ls)-3
    checks_passed = input('There are less saved data than contained in the Bacmman folder.\nAre you sure you wish to continue and erase all current Bacmman data for this folder?\n(y/n)\n', 's');
end
if checks_passed ~= 'y'
    return
end

% Perform saving and/or loading of data
for i = 1:length(ls)
    if ~sum(strfind('..Selections', ls(i).name))
        
        clc
        disp(['Processing ' num2str(i) '/' num2str(length(ls))])
        
        % Save folders
        if save_preprocessing
            copyfile([bacmman_folder '/Output/' ls(i).name '/pre_processed'], [save_folder '/Output/' ls(i).name '/pre_processed'])
        end
        if save_segmentation
            copyfile([bacmman_folder '/Output/' ls(i).name '/segmented_objects'], [save_folder '/Output/' ls(i).name '/segmented_objects'])
        end
        
        % Load folders
        if load_preprocessing
           rmdir([bacmman_folder '/Output/' ls(i).name '/pre_processed'], 's') % Remove previous pre-processing from Bacmman folder
           copyfile([save_folder '/Output/' ls(i).name '/pre_processed'], [bacmman_folder '/Output/' ls(i).name '/pre_processed']) % Copy saved pre-processing
        end
        if load_segmentation
            rmdir([bacmman_folder '/Output/' ls(i).name '/segmented_objects'], 's') % Remove previous segmentation from Bacmman folder
            copyfile([save_folder '/Output/' ls(i).name '/segmented_objects'], [bacmman_folder '/Output/' ls(i).name '/segmented_objects']) % Copy saved segmentation
        end
        
    end
end

clc
disp('Done!')


