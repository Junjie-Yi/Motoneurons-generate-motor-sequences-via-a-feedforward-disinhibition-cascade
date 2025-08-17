% This script allows the user to select a folder of images,
% and for each image, draw two lines and calculate the angle
% between them. The results are saved to an Excel file.

% 1. Ask the user to select a folder containing the images.
% 'uigetdir' opens a dialog to choose a directory.
% 'pwd' provides the current directory as a starting point.
folder = uigetdir(pwd, 'Select the folder containing images');

% Check if the user cancelled the selection.
if folder == 0
    disp('Folder selection cancelled. The script will now exit.');
    return;
end

% 2. Get a list of all image files in the selected folder.
% We assume '.bmp' files, but this can be changed to '*.jpg', '*.png', etc.
% 'fullfile' is used to construct the path, making it OS-independent.
image_files = dir(fullfile(folder, '*.bmp'));

% Check if any image files were found.
if isempty(image_files)
    disp('No .bmp images found in the selected folder. The script will now exit.');
    return;
end

% 3. Initialize a table to store the results.
% A table is a good way to organize data with named columns.
num_images = length(image_files);
results = table('Size', [num_images, 2], ...
                'VariableTypes', {'string', 'double'}, ...
                'VariableNames', {'FileName', 'Angle'});

% 4. Loop through each image file found in the folder.
for i = 1:num_images
    % Get the full path of the current image.
    file_name = image_files(i).name;
    full_file_path = fullfile(folder, file_name);

    % Read and display the image.
    img = imread(full_file_path);
    fig = figure('Name', ['Draw on: ' file_name], 'NumberTitle', 'off');
    imshow(img);

    % 5. Allow the user to select two points to define each line.
    disp(['Processing image: ' file_name]);
    title({'Select two points for the first line.'});
    disp('Please select two points for the first line.');
    [x1, y1] = ginput(2);
    
    % Draw the line for visual feedback
    hold on;
    line(x1, y1, 'Color', 'red', 'LineWidth', 2);
    hold off;
    
    title({'Select two points for the second line.'});
    disp('Please select two points for the second line.');
    [x2, y2] = ginput(2);

    % Draw the second line for visual feedback
    hold on;
    line(x2, y2, 'Color', 'green', 'LineWidth', 2);
    hold off;

    % Wait for the user to confirm.
    title({'Press any key to calculate angle and proceed.'});
    disp('Press any key in the figure to continue...');
    pause;

    % 6. Get the coordinates of the endpoints of the two lines.
    pos1 = [x1, y1];
    pos2 = [x2, y2];
    
    % Close the figure for the current image.
    close(fig);

    % 7. Calculate the angle between the two lines.
    % Create vectors from the line coordinates.
    v1 = [pos1(2,1) - pos1(1,1), pos1(2,2) - pos1(1,2)];
    v2 = [pos2(2,1) - pos2(1,1), pos2(2,2) - pos2(1,2)];

    % Calculate the signed angle from v1 to v2.
    % We use atan2 on the 2D cross product and the dot product.
    % atan2(v1(1)*v2(2)-v1(2)*v2(1), dot(v1,v2)) gives the angle in radians.
    % The result is positive for counter-clockwise and negative for clockwise.
    signed_angle_rad = atan2(v1(1)*v2(2) - v1(2)*v2(1), dot(v1, v2));

    % The user wants clockwise angles to be positive and counter-clockwise negative.
    % So we negate the result of atan2.
    angle_deg = -rad2deg(signed_angle_rad);

    % Display the calculated angle.
    fprintf('The calculated angle for %s is: %.2f degrees\n', file_name, angle_deg);

    % 8. Store the file name and the calculated angle in the results table.
    results.FileName(i) = string(file_name);
    results.Angle(i) = angle_deg;
end

% 9. Save the results to an Excel file.
% The file will be named based on the folder that was processed.
[~, folder_name, ~] = fileparts(folder);
output_filename = fullfile(folder, ['angle_results_' folder_name '.xlsx']);
try
    writetable(results, output_filename);
    fprintf('Results successfully saved to: %s\n', output_filename);
catch ME
    fprintf('Could not save the results to Excel file.\nError: %s\n', ME.message);
    % If saving to excel fails, save to a .mat file as a backup
    output_filename_mat = fullfile(folder, ['angle_results_' folder_name '.mat']);
    save(output_filename_mat, 'results');
    fprintf('Results saved to a MAT file instead: %s\n', output_filename_mat);
end

disp('All images have been processed.');
