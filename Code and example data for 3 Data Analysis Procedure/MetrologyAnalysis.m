clear; close all; clc;

intercept = true;
%% Read info
info = xlsread('info.xlsx');
dim = info(1); % Dimension along which the data was measured (1 for X, 2 for Y)
h = info(2); % Layer height
n = info(3); % Number of layers per step
Nsteps = info(4); % Number of vertical steps in the structure
edge_start = info(5); % Number of data points to be discarded at the start of the dataset (edge of the sample)
edge_end = info(6); % Number of data points to be discarded at the end of the dataset
stepMargin = info(7); % Number of data samples that are considered part of a step transition, to each side of the step
Zfraction = info(8); % Steps are found by searching for steep height increments of value [h*n*(1-1/Zfraction)] within a short data span (defined by 2*stepMargin)
Nz = info(9); % Number of data samples taken per file, per step, for computation of averages and standard deviations of step heights
Nxy = info(10); % Number of data samples taken per file, per step, for computation of averages and standard deviations of step widths
stepIncr = info(11); % In-plane step increment, as defined in CAD (used for comparison)
topWidth = info(12); % Width of top step, as defined in CAD (used for comparison)
smoothener = info(13); % Number of data points used to smoothen the dataset in order to identify steps
csvType = info(14); % Type of csv file format

%% Find files
csvFiles = dir('*.csv');
data = struct([]);
% Create results variables
sampled_h = zeros(Nz*length(csvFiles),2*Nsteps);
pos_sampled_h = zeros(Nz*length(csvFiles),2*Nsteps);
sampled_xy = zeros(Nxy*length(csvFiles),Nsteps);

%% Read files and prepare data structure
figure
hold
grid on
grid minor
for i = 1:length(csvFiles)
    data(i).name = csvFiles(i).name;
    switch csvType
        case 1
            data_aux = fileread(csvFiles(i).name); % Read file
            data_aux = strrep(data_aux,',','.'); % Change commas to dots
            newfile = fopen('dots.csv','w'); % Write new file with dots
            fwrite(newfile,data_aux,'char');
            fclose(newfile);
            data(i).profile = table2array(readtable('dots.csv','VariableNamingRule','preserve')); % Read new file (with dots)
            delete('dots.csv');
        case 2
            data(i).profile = readmatrix(csvFiles(i).name); % Read file
        otherwise
            disp('Unknown csv file format');
    end
    % Plot raw data
    plot(data(i).profile(:,1),data(i).profile(:,2));
end

%% Analysis
steps_sum = zeros(1,length(2*Nsteps));
for i = 1:length(data)
    % Discard measurements around the edges
    data(i).profile = data(i).profile(edge_start+1:end-edge_end,:);
    % Search for steps
    % Determine in-plane (X/Y) limits of steps
    steps = zeros(1,2*Nsteps);
    steps_index = steps;
    % Smooth data to get rid of noise and outlayers
    data_smooth = smoothdata(data(i).profile(:,2),'movmedian',smoothener);
    j = 1;
    left = 1;
    right = 1;
    % Swipe along the dataset searching for steep transitions
    while j < length(data(i).profile(:,1))-2*stepMargin
        if data_smooth(j+2*stepMargin)-data_smooth(j) > h*n*(1-1/Zfraction)
            % Rising edges
            steps(left) = data(i).profile(j+2*stepMargin,1);
            steps_index(left) = j+2*stepMargin;
            left = left + 1;
            j = j + 4*stepMargin;
        elseif data_smooth(j+2*stepMargin)-data_smooth(j) < -h*n*(1-1/Zfraction)
            % Falling edges
            steps(Nsteps+right) = data(i).profile(j+2*stepMargin,1);
            steps_index(Nsteps+right) = j+2*stepMargin;
            right = right + 1;
            j = j + 4*stepMargin;
        else
            j = j + 1;
        end
    end
    % Divide top step in two, so it weights as much as the other heights
    mid_step_index = round(mean([steps_index(Nsteps) steps_index(Nsteps+1)]));
    steps_index = [steps_index(1:Nsteps) mid_step_index steps_index(Nsteps+1:end)];
    steps = [steps(1:Nsteps) data(i).profile(mid_step_index,1) steps(Nsteps+1:end)];
    % Check that all steps have been identified, show warning if not
    if left < Nsteps+1
        disp(strcat('Not enough risign steps were found on dataset ',data(i).name));
    end
    if left > Nsteps+1
        disp(strcat('Too many risign steps were found on dataset ',data(i).name));
    end
    if right < Nsteps+1
        disp(strcat('Not enough falling steps were found on dataset ',data(i).name));
    end
    if right > Nsteps+1
        disp(strcat('Too many falling steps were found on dataset ',data(i).name));
    end
    data(i).steps = steps;
    % Reference heights to average bottom height
    bottom = mean([data(i).profile(1:(steps_index(1)-2*stepMargin),2); data(i).profile((steps_index(end)+2*stepMargin):end,2)]);
    data(i).profile(:,2) = data(i).profile(:,2) - bottom;
    data_smooth = data_smooth - bottom;

    % Determine out-of-plane (Z) limits of steps
    % Using the step edges in X/Y computed previously, the Z data within
    % edges is averaged to find approximate step transitions in Z
    zsteps = zeros(1,2*Nsteps);
    for j = 1:2*Nsteps
        zsteps(j) = mean(data_smooth(find(data(i).profile(:,1)>steps(j) & data(i).profile(:,1)<steps(j+1))));
    end
    
    % Plot data and identified step edges
    figure
    hold
    plot(data(i).profile(:,1),data(i).profile(:,2))
    plot(data(i).profile(:,1),data_smooth)
    xline(steps,'r')
    yline(zsteps,'g')
    
    %% Data extraction
    % Out-of-plane (Z) sampling
    % Extract Nz data samples for the height of each step (samples of Z)
    % A margin of 'stepMargin' samples on each side of the steps is discarded as part
    % of the step transition. Nz samples within the remaining data are
    % collected
    for j = 1:length(steps)-1
        index = round(linspace(steps_index(j)+stepMargin,steps_index(j+1)-stepMargin,Nz));
        sampled_h(Nz*(i-1)+1:Nz*i,j) = data(i).profile(index,2);
        pos_sampled_h(Nz*(i-1)+1:Nz*i,j) = data(i).profile(index,1);
        % If a chosen datapoint is clearly an outlayer, it is substituted
        % by a sample taken from the smoothed dataset
        for k = 1:Nz
            if sampled_h(Nz*(i-1)+k,j) < (data_smooth(index(k)) - 1.2*h*n) || sampled_h(Nz*(i-1)+k,j) > (data_smooth(index(k)) + 1.2*h*n)
                sampled_h(Nz*(i-1)+k,j) = data_smooth(index(k));
            end
        end
    end
    plot(pos_sampled_h(Nz*(i-1)+1:Nz*i,:),sampled_h(Nz*(i-1)+1:Nz*i,:),'*')
    
    % In-plane (X or Y) sampling
    for j = 1:Nsteps
        % Define valid sampling height range for each step
        if j == 1
            range_low = 0;
        else
            range_low = max([zsteps(j-1) zsteps(2*Nsteps+2-j)]);
        end
        range_high = min([zsteps(j) zsteps(2*Nsteps+1-j)]);
        % Choose Nxy sampling heights within the range, with some margin
        sample_z = linspace(range_low+(range_high-range_low)/5,range_high-(range_high-range_low)/5,Nxy);
        yline(sample_z,':c')
        % On each side of the step, find the data points closest to the
        % selected sampling heights, within a certain margin from the step
        % edge in XY
        % Left side
        A = repmat(data(i).profile(steps_index(j)-stepMargin:steps_index(j)+stepMargin,2),[1 length(sample_z')]);
        [minValueA, closestIndexA] = min(abs(A-sample_z));
        closestIndexA = steps_index(j)-stepMargin+closestIndexA-1;
        closestValueA = data(i).profile(closestIndexA,2);
        plot(data(i).profile(closestIndexA,1),data(i).profile(closestIndexA,2),'*c')
        % Right side
        B = repmat(data(i).profile(steps_index(2*Nsteps+2-j)-stepMargin:steps_index(2*Nsteps+2-j)+stepMargin,2),[1 length(sample_z')]);
        [minValueB, closestIndexB] = min(abs(B-sample_z));
        closestIndexB = steps_index(2*Nsteps+2-j)-stepMargin+closestIndexB-1;
        closestValueB = data(i).profile(closestIndexB,2);
        plot(data(i).profile(closestIndexB,1),data(i).profile(closestIndexB,2),'*c')
        % Compute in-plane difference
        XYdiff = data(i).profile(closestIndexB,1) - data(i).profile(closestIndexA,1);
        sampled_xy(Nxy*(i-1)+1:Nxy*i,Nsteps+1-j) = XYdiff;
    end

end

% Calculate averages and deviations
% Out-of-plane
avg_h = mean(sampled_h);
std_h = std(sampled_h);
figure
errorbar(avg_h,std_h,'o')
xlim([0 2*(Nsteps)+1])
yticks(0:h*n:h*n*Nsteps)
xlabel('Steps')
ylabel('Height (um)')
title('Height averages and standard deviations')
grid on
% In-plane
avg_xy = mean(sampled_xy);
std_xy = std(sampled_xy);
figure
errorbar(avg_xy,std_xy,'o')
CAD_xy = topWidth:2*stepIncr:topWidth+(Nsteps-1)*2*stepIncr;
xlim([0 Nsteps+1])
yticks(CAD_xy)
xlabel('Steps from the top')
ylabel('Width (um)')
title('Step width averages and standard deviations')
grid on

%% CAD vs measurements comparison
% Out-of-plane
heights = [sampled_h(:,1:Nsteps); flip(sampled_h(:,Nsteps+1:end),2)];
CADheights = ones(length(heights),1)*(n*h:n*h:n*h*Nsteps);
heights_line = reshape(heights,[],1);
CADheights_line = reshape(CADheights,[],1);
lineZ = fitlm(CADheights_line,heights_line);
figure
plot(lineZ)
xlim([0 n*h*(Nsteps+1)])
yticks(0:h*n:h*n*Nsteps)
grid on
title('Straight line fit to measured heights vs CAD heights')
xlabel('CAD dimensions (um)')
ylabel('Measured dimensions (um)')
zline = strcat('z=',num2str(lineZ.Coefficients.Estimate(1)),'+',num2str(lineZ.Coefficients.Estimate(2)),'*zCAD');
zRsq = strcat('R^2=',num2str(lineZ.Rsquared.Adjusted));
zmodel = {zline,zRsq};
text(n*h,n*h*(Nsteps-2),zmodel)

% In-plane
switch dim
    case 1
        dimension = 'x';
    case 2
        dimension = 'y';
end
CADwidths = ones(length(sampled_xy),1)*CAD_xy;
widths_line = reshape(sampled_xy,[],1);
CADwidths_line = reshape(CADwidths,[],1);
lineXY = fitlm(CADwidths_line,widths_line,'Intercept',intercept);
figure
plot(lineXY)
xlim([topWidth-2*stepIncr topWidth+Nsteps*2*stepIncr])
xticks(CAD_xy)
yticks(CAD_xy)
grid on
title('Straight line fit to measured step widths vs CAD widths')
xlabel('CAD dimensions (um)')
ylabel('Measured dimensions (um)')
if intercept
    xyline = strcat(dimension,'=',num2str(lineXY.Coefficients.Estimate(1)),'+',num2str(lineXY.Coefficients.Estimate(2)),'*',dimension,'CAD');
else
    xyline = strcat(dimension,'=',num2str(lineXY.Coefficients.Estimate(1)),'*',dimension,'CAD');
end
xyRsq = strcat('R^2=',num2str(lineXY.Rsquared.Adjusted));
xymodel = {xyline,xyRsq};
text(topWidth,topWidth+2*stepIncr*(Nsteps-2),xymodel)
% Taking out the top step
sampled_xy_cut = sampled_xy(:,2:end);
CAD_xy_cut = CAD_xy(:,2:end);
CADwidths_cut = ones(length(sampled_xy_cut),1)*CAD_xy_cut;
widths_line_cut = reshape(sampled_xy_cut,[],1);
CADwidths_line_cut = reshape(CADwidths_cut,[],1);
lineXYcut = fitlm(CADwidths_line_cut,widths_line_cut);
figure
plot(lineXYcut)
xlim([topWidth-2*stepIncr topWidth+Nsteps*2*stepIncr])
xticks(CAD_xy_cut)
yticks(CAD_xy_cut)
grid on
title('Straight line fit to measured step widths vs CAD widths WITHOUT TOP STEP')
xlabel('CAD dimensions (um)')
ylabel('Measured dimensions (um)')
xyline_cut = strcat(dimension,'=',num2str(lineXYcut.Coefficients.Estimate(1)),'+',num2str(lineXYcut.Coefficients.Estimate(2)),'*',dimension,'CAD');
xyRsq_cut = strcat('R^2=',num2str(lineXYcut.Rsquared.Adjusted));
xymodel_cut = {xyline_cut,xyRsq_cut};
text(topWidth,topWidth+2*stepIncr*(Nsteps-2),xymodel_cut)

%% Save data
results.heights = heights_line;
results.heightsCAD = CADheights_line;
results.widths = widths_line;
results.widthsCAD = CADwidths_line;
results.widths_cut = widths_line_cut;
results.widthsCAD_cut = CADwidths_line_cut;
results.dim = dimension;
results.info = strcat(num2str(h),'-',num2str(n),'-',dimension);
dir = pwd;
sep = strfind(dir,'\');
save(strcat(dir(sep(end-1)+1:sep(end)-1),dir(sep(end)+1:end)),'results');
