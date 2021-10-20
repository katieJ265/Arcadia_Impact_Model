%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculates change in labour productivity                            %
%  1. Calcualte # Days which exceed WGBT threshlolds for different     %
%       work intesities (with and without acclimatisation)             %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUTS:                                                             %
%   1. labour_productivity_ERFs [4,5]: sWGBT Thresholds below which    %
%      which labour productivity is unaffected (row 1); declines       %
%      between given lower and upper sWGBT range(rows 1-2) based on    %
%      slope of trend (a and b parameters, rows 3-4) and above which   %
%      100% decline in labour productivity (row 2). Assumes no         %
%      acclimatisation. Based on ERFs in Costa et al. (2016). Columns  %
%      provide values per 5 work intensity categories.                 %
%   2. labour_productivity_ERFs_acc [4,5]: Same as above but assume    %
%      acclimatised workers as default.                                %
%   3. workIntensity: [16,1] Defines the average work intensity per    %
%      sector from 1-5, 5 being high and 1 being light), expanded from %
%      Costa et al. (2016). Defines which ERF to use.
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp ('Running labour productivity calculations')

%% 1. Search for available files; Read temperature data from HEAT .csv files; set up arrays. 
%  Read each gridcell array to a cell Array inc. the cell ID.
%  Each TMean array represents daily data e.g. [1-360] per year [1-30][nCols=10,800] * 12 GCMS [nRows = 12]
%  Only need to do once - if files exist then skip

% RUN EITHER BIAS CORRECTED OR RAW DEPENDENT ON USER DEFINED PARAMETER

if BiasCorrected == 0 %False, reads raw data which has 30 day months, no leap years
   disp ('Using non-bias corrected data')
  
      %  Past sWGBT
      %  If raw data already read in once and saved then skip initial step
    if exist ('input_data_past_wgbt.mat', 'file')
       load ('input_data_past_wgbt.mat')
       disp ('loading input data Past sWGBT')
       nClimScen = nClimScen+1; %counts # climate scenarios being run together

    elseif exist('ARCADIA_land_RCM_sWGBT_raw_past', 'file') 
        disp ('Read in daily sWGBT data (Past) for labour productivity estimates')

        dd = dir('ARCADIA_land_RCM_sWGBT_raw_past\*.csv'); % all .csv files in folder
        fileNames = {dd.name}; 
        data_past_wgbt = cell(numel(fileNames),2);
        data_past_wgbt(:,1) = regexprep(fileNames, '_UKCP18-12km_sWBGT.csv',''); %set name as grid_cell ID (long_lat)

        for ii = 1:numel(fileNames)   
            fname = fullfile('ARCADIA_land_RCM_sWGBT_raw_past', fileNames(ii)); %allocate each .csv data file to cell array
            data_past_wgbt{ii,2} = dlmread(fname{1},',',0,0);
        end

        save('input_data_past_wgbt.mat', 'data_past_wgbt', '-v7.3') %~1.1GB (cell_IDs and data array)
        nClimScen = nClimScen+1; %counts # climate scenarios being run together

    else
        disp ('No Past files exist')
    end
end

%% Set lat and long for use in mapping outputs and linking regional mortality thresholds and Relative Risk (RR) to gridded TMean
%(format DS --> currently convert to DMS in ArcGIS for plotting, but could convert here if required as output in this format)
% Only need to do once, if file exists skip step.

if exist ('gridIDwgbt.mat', 'file')
    load ('gridIDwgbt')  
    disp('loading WGBT gridID data')
    [nGridCells,nVar] = size(gridIDwgbt);
else
    disp('calculating gridID data')
    load lat_UK_RCM
    load long_UK_RCM
    [nGridCells,nVar] = size(data_past_wgbt);
    gridIDwgbt = zeros(nGridCells,2); % sets corresponding lon (x) lat (y) for gridCells for plotting

    cellIndex = split(data_past_wgbt(:,1),"_"); % remove hyphen. cellIndex used to pull out corresponding land-based long_lat from UK & ROI grid
    cellIndex = str2double(cellIndex);

    for i = 1:nGridCells
        x = cellIndex(i,1);
        y = cellIndex(i,2);
        gridIDwgbt(i,2) = lat_UK_RCM(x,y); %corresponding long_lat for each cell read in. Corresponds to rows of [data_climScen] and [population].
        gridIDwgbt(i,1) = long_UK_RCM(x,y);
    end
    gridIDwgbt = fix(gridIDwgbt*1e5)/1e5;
   
    %UPDATE USING EMPLOYMENT GRIDDED DATA
    lonIDIndexWBGT = fix(Mortality_Threshold_RR(:,1)*1e5)/1e5; %truncate decimals to allow match to gridID - different order to IDs in gridID
    latIDIndexWBGT = fix(Mortality_Threshold_RR(:,2)*1e5)/1e5; %truncate decimals to allow match to gridID - different order to IDs in gridID
    
    save('gridIDwgbt', 'gridIDwgbt', 'lonIDIndexWBGT', 'latIDIndexWBGT')
    clear cellIndex %not used past this point - use gridID
end

%% Define size of variables / Cell Arrays for calculations
    [nRows,nCols] = size(data_past_wgbt{1,2}); %set up TMean data array size
    years = nCols/daysPerYear; %check years correspond to 30
    nSector = size(workIntensity, 1); %check number of sectors data read in for (default 16)
    nWorkIntensity = size(labour_productivity_ERFs,2); % cthe number of work intesnity categories data read in for (default 5)
    
%% Calculate percentage change in productiviyt
% For each day see if the value of TMean (will be WGBT) exceeds the
% thresholds and calculate annual aggregate reduction in productivity.

acclimatised = 1; % 0 = non-acclimatised; 1= acclimatised. Defines which Exposure response Functions (ERFs) are used.
if acclimatised == 1
    ERFs = labour_productivity_ERFs_acc;
elseif acclimatised == 0
    ERFs = labour_productivity_ERFs;
end

h=1; %1:nGridCells
nWorkIntensityGroups = 5; %WI1-5 acclimatised, WI6-10 non-acclimatised
dailyHrs = 8; % assumes typical 8 hour work day
AvProdLossHrs = zeros(nRows,nGridCells);
AvProdLossPercent = zeros(nRows,nGridCells);
a = 1; %%update if have different adaptation scenarios in future.
hrsPerYear = dailyHrs*264; %assume work 22 days per month on average, 264
%per year.
% nSectors = size(employees,1); 

%% PAST
%Check if results already available - only run once then save.
if BiasCorrected == 0
    if exist ('ResultsPastLabour.mat', 'file')
       load ('ResultsPastLabour')
       disp ('loading Results: Past for plotting')
       
    elseif exist ('input_data_past_wgbt.mat', 'file')  %calculate for each time period - based on input files %CHANGE to WGBT
    disp ('Calculating change in labour productivity - Near Past')
    
    
    %Preallocate cell arrays
    ResultsPastLabour = cell(adaptScen,1); %append all final gridded results for plotting and saving, for each age group, one array for each adaptation scenario
    for n = 1:adaptScen
        ResultsPastLabour{n} = zeros(nGridCells, 17);
    end

    ResultsPastLabourUK = cell(adaptScen,1); %append all final results aggregated per age group, for plotting and saving, one array for each adaptation scenario
    for n = 1:adaptScen
        ResultsPastLabourUK{n} = zeros(5, 3); %%%replace with nSector when read in employment data.
    end

    annProdLossHrs = cell(nGridCells,1);   %store intermediate value of change in hours loss - aggregated to annual loss.       
        for n = 1:nGridCells
            annProdLossHrs{n} = zeros(nRows, nCols);
        end
    
    % Replace gridded population with gridded # employees per sector
    %Calculate 
    for w = 1:nWorkIntensityGroups
        
        changeLabourProd = cell(nGridCells,1);   %store intermediate value of daily change in labour productivity (%)     
        for n = 1:nGridCells
            changeLabourProd {n} = zeros(nRows, nCols);
        end
           %calculates the loss in hours worked daily for different work
           %intensities
        
          i=1;
          while (h <= nGridCells)
               if gridIDwgbt(i,1) == lonIDIndexWBGT(h,1) && gridIDwgbt(i,2) == latIDIndexWBGT(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                  for r = 1:nRows                                                          % i.e. makes sure each lon/lat index, used in file name, is assigned the correct DS Grid ID coordinate for plotting.
                    for c = 1:nCols
                        if (data_past_wgbt{i,2}(r,c) < ERFs(1,w))
                            changeLabourProd{h}(r,c) = 0;
                        elseif(data_past_wgbt{i,2}(r,c) >= ERFs(1,w)) && (data_past_wgbt{i,2}(r,c) <= ERFs(2,w))    %less than value results remains zero
                            changeLabourProd{h}(r,c) = dailyHrs*(1-(ERFs(3,w)-ERFs(4,w)*data_past_wgbt{i,2}(r,c)));
                        elseif (data_past_wgbt{i,2}(r,c) > ERFs(2,w))
                            changeLabourProd{h}(r,c) = dailyHrs; %i.e. 100% loss in productivity on that day
                        end
                           
                    end
                  end
                    h=h+1; %increment through data in Mortality_Threshold_RR
                    i=1; %reset counter
                    r=1; %reset counter
                    c=1; %reset counter
               else
                 i=i+1;
               end
           end
           
           %calculates average loss in hours per year. Sumas daily data to year -> total hrs lost per year for 30 years
           for h = 1:nGridCells
                annProdLossHrs{h} = squeeze(nansum(reshape(changeLabourProd{h}.',daysPerYear,size(changeLabourProd{h},2)./daysPerYear,[]))).'; %re-shape array
                for r = 1:nRows
                    AvProdLossHrs(r,h) = mean(annProdLossHrs{h}(r,:),[1,2]); %mean across 30year period for each of 12 RCM (r) per cell (h)
                    %AvProdLossPercent(r,h) = (mean(annProdLossHrs{h}(r,:),[1,2]))/hrsPerYear*100; %average annual % loss in labour productivity for each of 12 RCM (r) per cell (h)
                end
                
                %summary values to save
                p_low = prctile(AvProdLossHrs,(10),1);%10th percentile across years/RCMs
                p_high = prctile(AvProdLossHrs,(90),1);%90th percentile across years/RCMs
                av = mean(AvProdLossHrs,1);%average across 12 RCMs
                sum_av = sum(av);
                sum_p_low = sum(p_low);
                sum_p_high = sum(p_high);
           end
           
            %Save results for plotting - gridded
                    if w == 1 % all ages
                        ResultsPastLabour{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                        ResultsPastLabour{a}(1:h,1) = lonIDIndexWBGT(:,1); %longitude
                        ResultsPastLabour{a}(1:h,2) = latIDIndexWBGT(:,1); %latitude
                    %Save results for plotting - UK aggregate
                        ResultsPastLabourUK{a}(w, 1:3) = [sum_av.', sum_p_low.', sum_p_high.'];

                    elseif w==2
                        ResultsPastLabour{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                        ResultsPastLabourUK{a}(w, 1:3) = [sum_av.', sum_p_low.', sum_p_high.'];
                    elseif w==3
                        ResultsPastLabour{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                        ResultsPastLabourUK{a}(w, 1:3) = [sum_av.', sum_p_low.', sum_p_high.'];
                    elseif w==4
                        ResultsPastLabour{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                        ResultsPastLabourUK{a}(w, 1:3) = [sum_av.', sum_p_low.', sum_p_high.'];
                    elseif w==5
                        ResultsPastLabour{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                        ResultsPastLabourUK{a}(w, 1:3) = [sum_av.', sum_p_low.', sum_p_high.'];
                    end

 
                    %Re-set counters
                    h=1; %1:nGridCells
     end
                    % Files to SAVE
                    save('ResultsPastLabour.mat', 'ResultsPastLabour', 'ResultsPastLabourUK') % loss in hours per year (gridded) and agrgeated annual hours lost (UK-wide), averaged across 12 RCMs to give mean
    end
end

%% Plots
% Create some initial output figures

%1. Average annual hours of labour productivity lost per WI category
% (intermediate output) - Test using near-past data

figure;
x=1:w;
subplot(1,2,1);
% Create bar and set individual colours
b = bar(ResultsPastLabourUK{1}(:,1));

hold on
p_errhigh = ResultsPastLabourUK{1}(:,3)- ResultsPastLabourUK{1}(:,1); %90th p - average
p_errlow = ResultsPastLabourUK{1}(:,1)-ResultsPastLabourUK{1}(:,2); % average - 10th p)

er = errorbar(x, ResultsPastLabourUK{1}(:,1), p_errlow, p_errhigh, 'LineWidth',1);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

% Create ylabel
ylabel({'Heat related impact on labour productivity:', 'Average Annual Hours lost'});

% Create xlabel
xlabel('Work Intensity Category');


% Set the remaining axes properties
set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'1','2','3','4', '5'},'YGrid',...
    'on');

subplot(1,2,2);
% Create bar  - percentage loss of annual hours
ResultsPastLabourUKPercent = (ResultsPastLabourUK{1}(:,:)/hrsPerYear*100);
b = bar(ResultsPastLabourUKPercent(:,1));

hold on
p_errhigh = ResultsPastLabourUKPercent(:,3)- ResultsPastLabourUKPercent(:,1); %90th p - average
p_errlow = ResultsPastLabourUKPercent(:,1)-ResultsPastLabourUKPercent(:,2); % average - 10th p)

er = errorbar(x, ResultsPastLabourUKPercent(:,1), p_errlow, p_errhigh, 'LineWidth',1);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

% Create ylabel
ylabel({'Heat related impact on labour productivity:', 'Percent Average Annual Hours lost'});

% Create xlabel
xlabel('Work Intensity Category');

% Set the remaining axes properties
set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'1','2','3','4', '5'},'YGrid',...
    'on');

% Create title
suptitle('Climate Change Only - Near Past');