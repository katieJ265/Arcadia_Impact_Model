%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ARCADIA_HEAT_Impacts_model
%
%   Estimates heat related mortality risk, residential discomfort, and
%   labour productivity risk.
%   
%   Updated from model developed in ARCADIA project (ECI, University of Oxford)
%   for OpenCLIM project.
%
%   Uses UKCP18 data provided from HEAT model (University of Bristol)
%
%   Author: Katie Jenkins    UEA
%   Date: 12/11/2020 [version_01]
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT DATA
%   heat_impact_data.mat
%   1.  Mortality_Threshold_RR (nGridCells * 9)
%       Column 1 and 2: Lon and Lat of grid cell
%       Column 3: NUTS1/GOR name index 1-13 (inc. RoI) (see GOR_names)
%       Column 4: Threshold TMean above which heat related mortality will
%       occur by GOR. Source: Vardoulakis et al., (2014) and Hajat et al., (2014)
%       Column 5-9: Exosure-Risk relationship (% change in mortality for
%       each 1degeeC increase in Mean Temp above the MortalityThreshold for
%       all ages and by age group (see Age_group). Source: Vardoulakis et al., 
%       (2014) and Hajat et al., (2014)).
%       
%   2.  Gridded population data {nClimScen}[nGridCells, 7]
%       Column 1&2 lat and long
%       column 3 {1} = Baseline 2011 census based gridded data (1km) 
%       regridded to 12km in ArcGIS. Source: Reis, S. et al.(2017) 
%       https://doi.org/10.5285/0995e94d-6d42-40c1-8ed4-5090d82471e1
%       column 4-7 {1} = gridded poulation per age group based on ONS
%       demographic projections from 2011 census for LAs.
%       Columns 4-7 {2-4} = For 2020, 2030, 2050 UK population projections and demographic data from
%       UK-SSPs SSP5. All ages 1km re-gridded to 12km. Per age group LAs gridded to 12km 
%       Source: https://www.ukclimateresilience.org/products-of-the-uk-ssps-project/
%       
%   3.  population_aggregate [5, 4] - for each climate scenario (columns)
%       gives the aggegrate UK population and split by age groups (rows).
%       Sources as above.
%
%   4.  Gridded dailyDeathRate calculated for 12km grid based on
%       mortality statistics from 2011 for England, Wales, Scotland and NI.
%       
%   5.  GOR_names [13,1]: North East, North West, Yorkshire and Humber, East Midlands,
%       West Midlands, East England, London, South East, South West, Wales,
%       Scotland, Northern Ireland, Republic of Ireland.
%
%   6.  Age_group [5,1]: description of the age group classes used.
%      ['All';'0-64';'65-74';'75-84';'85+']
%
%   7. 12km gridded TMean daily data {nGridCells}(12, 14800). .csv files of 
%       daily TMean (one per grid cell) from HEAT model for past (1990-2019); and warming levels 1.5°C (1995-2024)
%       2.0°C(2005-2034) and 3.0°C (2023-2052) 
%
%       These are read in, saved and used for mortality calculations. 
%       Grid coordinates also extracted from file names to provide grid ID for 
%       mapping based on lat_UK_RCM.mat and long_UK_RCM.mat.
%
%   8/9. lat_UK_RCM.mat and long_UK_RCM.mat: [82,112]. From Heat outputs, these are read 
%      in and provide latitude and longitude based on coordinate given in file 
%      names from 6 above. Saved in gridID.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CURRENT MODEL OUTPUTS:
%   1. Plots average annual heat related deaths (10th & 90th P reflect climate
%   model uncertainty), for CC (Climate change Only Scenario) and CC+SE
%   (Climate Change and future projection of population both accounted for).
%   Scenarios reflect the past, 1.5, 2.0 and 3 degree warming levels.
%   2. Heat related deaths per year per 100,000 population (Total and split
%   by age groups)
%   3. Annual het related deaths, all years*12 RCM - for looking at whole
%   distribution rather than just average annual values.
%   4. Daily and annual timeseries as intermediate outputs if needed, split by GCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('*******************************************************************')
disp('ARCADIA Heat related Mortality Impacts Model')
disp('version 01')
disp('This version currently only estimates heat-related mortality')
disp('*******************************************************************')

tic %stopwatch to measure performance
format short e
dbstop if error %debug if error

load('heat_impact_data') %contains all the .m input data required.
disp ('loading input parameters')

%Parameters
daysPerYear = 0; %UKCP18 has 360 days per year (30 days per month)- set to 360 for raw data below %Bias corrected converted to 365 days per year (leap years removed below) - set to 365 for bias corrected data below
nClimScen = 0; %number of climate scenarios being read in...
socioEcScen = 2; % Climate change only && climate change and socio-economic change
adaptScen = 3; %Used to calculate results for acclimatisation scenarios by adjusting MortalityThreshold: None = 1; 1degreeC = 2; 2degreeC =3
cnt = 1; % counter to loop through acclimatisation adaptation scenarios(adaptScen), i.e. +1 and +2degreeC acclimatisation scenario

%User defined parameters - change here
BiasCorrected = 0; %TRUE = 1, FALSE = 0 %%USER DEFINED AT THIS STAGE
adaptIncrement = [0,1,2]; %i.e. +1 and +2degreeC acclimatisation scenario
impactMetric = 2; %Select impact to run - 1 = mortality; 2 = labour productivity

if BiasCorrected == 0 %false
        daysPerYear = 360;
elseif BiasCorrected == 1 %true
        daysPerYear = 365;
end

if impactMetric ==1  %%Mortality calculations
    
%% 1. Search for available files; Read temperature data from HEAT .csv files; set up arrays. 
%  Read each gridcell array to a cell Array inc. the cell ID.
%  Each TMean array represents daily data e.g. [1-360] per year [1-30][nCols=10,800] * 12 GCMS [nRows = 12]
%  Only need to do once - if files exist then skip

% RUN EITHER BIAS CORRECTED OR RAW DEPENDENT ON USER DEFINED PARAMETER

if BiasCorrected == 0 %False, reads raw data which has 30 day months, no leap years
   disp ('Using non-bias corrected data')
 
    %  Past (1990-2019) TMean
    %  If raw data already read in once and saved then skip initial step
    if exist ('input_data_past.mat', 'file')
       load ('input_data_past.mat')
       disp ('loading input data Past TMean')
       nClimScen = nClimScen+1; %counts # climate scenarios being run together

    elseif exist('ARCADIA_land_RCM_raw_past', 'file') 
        disp ('Read in daily TMean data (Past) for mortality estimates')

        dd = dir('ARCADIA_land_RCM_raw_past\*.csv'); % all .csv files in folder
        fileNames = {dd.name}; 
        data_past = cell(numel(fileNames),2);
        data_past(:,1) = regexprep(fileNames, '_UKCP18-raw_past.csv',''); %set name as grid_cell ID (long_lat)

        for ii = 1:numel(fileNames)   
            fname = fullfile('ARCADIA_land_RCM_raw_past', fileNames(ii)); %allocate each .csv data file to cell array
            data_past{ii,2} = dlmread(fname{1},',',0,0);
        end

        save('input_data_past.mat', 'data_past', '-v7.3') %~1.1GB (cell_IDs and data array)
        nClimScen = nClimScen+1; %counts # climate scenarios being run together

    else
        disp ('No Past files exist')
    end

    %  1.5degreeC (1995-2024)
    if exist ('input_data_1.5C.mat', 'file')
       load ('input_data_1.5C.mat')
       disp ('loading input data 1.5C')
       nClimScen = nClimScen+1; %counts # climate scenarios being run together

    elseif exist('ARCADIA_land_RCM_raw_1_5C', 'file')
        disp ('Read in daily TMean data (1.5DegreeC) for mortality estimates')

        dd = dir('ARCADIA_land_RCM_raw_1_5C\*.csv'); % all .csv files in folder
        fileNames = {dd.name}; 
        data_1_5C = cell(numel(fileNames),2);
        data_1_5C(:,1) = regexprep(fileNames, '_UKCP18-raw_1.5C.csv',''); %set name as grid_cell ID (long_lat)

        for ii = 1:numel(fileNames)
            fname = fullfile('ARCADIA_land_RCM_raw_1_5C', fileNames(ii)); %allocate each .csv data file to cell array
            data_1_5C{ii,2} = dlmread(fname{1},',',0,0);
        end

        save('input_data_1.5C.mat', 'data_1_5C', '-v7.3') %~1.1GB (cell_IDs and data array)
        nClimScen = nClimScen+1; %counts # climate scenarios being run together

    else
        disp ('No 1.5degreeC files exist')
    end

    %  2degreeC (2005-2034)
    if exist ('input_data_2.0C.mat', 'file')
       load ('input_data_2.0C.mat')
       disp ('loading input data 2.0C')
       nClimScen = nClimScen+1; %counts # climate scenarios being run together

    elseif exist('ARCADIA_land_RCM_raw_2_0C', 'file')
        disp ('Read in daily TMean data (2.0DegreeC) for mortality estimates')

        dd = dir('ARCADIA_land_RCM_raw_2_0C\*.csv'); % all .csv files in folder
        fileNames = {dd.name}; 
        data_2_0C = cell(numel(fileNames),2);
        data_2_0C(:,1) = regexprep(fileNames, '_UKCP18-raw_2.0C.csv',''); %set name as grid_cell ID (long_lat)

        for ii = 1:numel(fileNames)   
            fname = fullfile('ARCADIA_land_RCM_raw_2_0C', fileNames(ii)); %allocate each .csv data file to cell array
            data_2_0C{ii,2} = dlmread(fname{1},',',0,0);
        end

        save('input_data_2.0C.mat', 'data_2_0C', '-v7.3') %~1.1GB (cell_IDs and data array)
        nClimScen = nClimScen+1; %counts # climate scenarios being run together

    else
        disp ('No 2.0DegreeC files exist')
    end

    %  3degreeC (2023-2051)
    if exist ('input_data_3.0C.mat', 'file')
       load ('input_data_3.0C.mat')
       disp ('loading input data 3.0C')
       nClimScen = nClimScen+1; %counts # climate scenarios being run together

    elseif exist('ARCADIA_land_RCM_raw_3_0C', 'file')
        disp ('Read in daily TMean data (3.0DegreeC) for mortality estimates')

        dd = dir('ARCADIA_land_RCM_raw_3_0C\*.csv'); % all .csv files in folder
        fileNames = {dd.name}; 
        data_3_0C = cell(numel(fileNames),2);
        data_3_0C(:,1) = regexprep(fileNames, '_UKCP18-raw_3.0C.csv',''); %set name as grid_cell ID (long_lat)

        for ii = 1:numel(fileNames)   
            fname = fullfile('ARCADIA_land_RCM_raw_3_0C', fileNames(ii)); %allocate each .csv data file to cell array
            data_3_0C{ii,2} = dlmread(fname{1},',',0,0);
        end

        save('input_data_3.0C.mat', 'data_3_0C', '-v7.3') %~1.1GB (cell_IDs and data array)
        nClimScen = nClimScen+1; %counts # climate scenarios being run together

    else
        disp ('No 3.0DegreeC files exist')
    end
    
elseif BiasCorrected == 1 %reads bias corrected data and removes leap years.
       disp ('Using bias corrected data')
    
    %  Past (1990-2019)
    if exist ('input_data_past_bias.mat', 'file')
       load ('input_data_past_bias.mat')
       disp ('loading input data Past (Bias corrected)')
       nClimScen = nClimScen+1; %counts # climate scenarios being run together

    elseif exist('ARCADIA_land_RCM_raw_past', 'file')
        disp ('Read in daily TMean data (Past) for mortality estimates')

        dd = dir('ARCADIA_land_RCM_bias_past\*.csv'); % all .csv files in folder
        fileNames = {dd.name}; 
        data_past = cell(numel(fileNames),2);
        data_past(:,1) = regexprep(fileNames, '_UKCP18-biascorr_past.csv',''); %set name as grid_cell ID (long_lat)

        for ii = 1:numel(fileNames)   
            fname = fullfile('ARCADIA_land_RCM_bias_past', fileNames(ii)); %allocate each .csv data file to cell array
            data_past{ii,2} = dlmread(fname{1},',',0,0);
            data_past{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
        end

        save('input_data_past_bias.mat', 'data_past', '-v7.3') %~1.1GB (cell_IDs and data array)
        nClimScen = nClimScen+1; %counts # climate scenarios being run together

    else
        disp ('No Past (bias corrected) files exist')
    end

    %  1.5degreeC (1995-2024)
    if exist ('input_data_1.5C_bias.mat', 'file')
       load ('input_data_1.5C_bias.mat')
       disp ('loading input data 1.5C (Bias corrected)')
       nClimScen = nClimScen+1; %counts # climate scenarios being run together

    elseif exist('ARCADIA_land_RCM_bias_1_5C', 'file')
        disp ('Read in daily TMean data (1.5DegreeC) for mortality estimates')

        dd = dir('ARCADIA_land_RCM_bias_1_5C\*.csv'); % all .csv files in folder
        fileNames = {dd.name}; 
        data_1_5C = cell(numel(fileNames),2);
        data_1_5C(:,1) = regexprep(fileNames, '_UKCP18-biascorr_1.5C.csv',''); %set name as grid_cell ID (long_lat)

        for ii = 1:numel(fileNames)
            fname = fullfile('ARCADIA_land_RCM_bias_1_5C', fileNames(ii)); %allocate each .csv data file to cell array
            data_1_5C{ii,2} = dlmread(fname{1},',',0,0);
            data_1_5C{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
        end

        save('input_data_1.5C_bias.mat', 'data_1_5C', '-v7.3') %~1.1GB (cell_IDs and data array)
        nClimScen = nClimScen+1; %counts # climate scenarios being run together

    else
        disp ('No 1.5degreeC files exist')
    end

    %  2degreeC (2005-2034)
    if exist ('input_data_2.0C_bias.mat', 'file')
       load ('input_data_2.0C_bias.mat')
       disp ('loading input data 2.0C (Bias corrected)')
       nClimScen = nClimScen+1; %counts # climate scenarios being run together

    elseif exist('ARCADIA_land_RCM_bias_2_0C', 'file')
        disp ('Read in daily TMean data (2.0DegreeC) for mortality estimates')

        dd = dir('ARCADIA_land_RCM_bias_2_0C\*.csv'); % all .csv files in folder
        fileNames = {dd.name}; 
        data_2_0C = cell(numel(fileNames),2);
        data_2_0C(:,1) = regexprep(fileNames, '_UKCP18-biascorr_2.0C.csv',''); %set name as grid_cell ID (long_lat)

        for ii = 1:numel(fileNames)   
            fname = fullfile('ARCADIA_land_RCM_bias_2_0C', fileNames(ii)); %allocate each .csv data file to cell array
            data_2_0C{ii,2} = dlmread(fname{1},',',0,0);
            data_2_0C{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
        end

        save('input_data_2.0C_bias.mat', 'data_2_0C', '-v7.3') %~1.1GB (cell_IDs and data array)
        nClimScen = nClimScen+1; %counts # climate scenarios being run together

    else
        disp ('No 2.0DegreeC files exist')
    end

    %  3degreeC (2023-2051)
    if exist ('input_data_3.0C_bias.mat', 'file')
       load ('input_data_3.0C_bias.mat')
       disp ('loading input data 3.0C (Bias corrected)')
       nClimScen = nClimScen+1; %counts # climate scenarios being run together

    elseif exist('ARCADIA_land_RCM_bias_3_0C', 'file')
        disp ('Read in daily TMean data (3.0DegreeC) for mortality estimates')

        dd = dir('ARCADIA_land_RCM_bias_3_0C\*.csv'); % all .csv files in folder
        fileNames = {dd.name}; 
        data_3_0C = cell(numel(fileNames),2);
        data_3_0C(:,1) = regexprep(fileNames, '_UKCP18-biascorr_3.0C.csv',''); %set name as grid_cell ID (long_lat)

        for ii = 1:numel(fileNames)   
            fname = fullfile('ARCADIA_land_RCM_bias_3_0C', fileNames(ii)); %allocate each .csv data file to cell array
            data_3_0C{ii,2} = dlmread(fname{1},',',0,0);
            data_3_0C{ii,2}(:,1461:1461:end) = []; %delete leap year, every 1461 columns. Bias corrected data uses Gregorian Calendar.
        end

        save('input_data_3.0C_bias.mat', 'data_3_0C', '-v7.3') %~1.1GB (cell_IDs and data array)
        nClimScen = nClimScen+1; %counts # climate scenarios being run together

    else
        disp ('No 3.0DegreeC files exist')
    end
       
end

%% Set lat and long for use in mapping outputs and linking regional mortality thresholds and Relative Risk (RR) to gridded TMean
%(format DS --> currently convert to DMS in ArcGIS for plotting, but could convert here if required as output in this format)
% Only need to do once, if file exists skip step.

if exist ('gridID.mat', 'file')
    load ('gridID')  
    disp('loading gridID data')
    [nGridCells,nVar] = size(gridID);
else
    disp('calculating gridID data')
    load lat_UK_RCM
    load long_UK_RCM
    [nGridCells,nVar] = size(data_past);
    gridID = zeros(nGridCells,2); % sets corresponding lon (x) lat (y) for gridCells for plotting

    cellIndex = split(data_past(:,1),"_"); % remove hyphen. cellIndex used to pull out corresponding land-based long_lat from UK & ROI grid
    cellIndex = str2double(cellIndex);

    for i = 1:nGridCells
        x = cellIndex(i,1);
        y = cellIndex(i,2);
        gridID(i,2) = lat_UK_RCM(x,y); %corresponding long_lat for each cell read in. Corresponds to rows of [data_climScen] and [population].
        gridID(i,1) = long_UK_RCM(x,y);
    end
    gridID = fix(gridID*1e5)/1e5;
   
    lonIDIndex = fix(Mortality_Threshold_RR(:,1)*1e5)/1e5; %truncate decimals to allow match to gridID - different order to IDs in gridID
    latIDIndex = fix(Mortality_Threshold_RR(:,2)*1e5)/1e5; %truncate decimals to allow match to gridID - different order to IDs in gridID
    
    save('gridID', 'gridID', 'lonIDIndex', 'latIDIndex')
    clear cellIndex %not used past this point - use gridID
end

%% Define size of variables / Cell Arrays for calculations
    [nRows,nCols] = size(data_past{1,2}); %set up TMean data array size
    years = nCols/daysPerYear; %check years correspond to 30
    ageGroupCnt = size(Age_group,1); %counter to loop through calculations for each age group
 
%% CALCULATE MORTALITY - Average Annual, per 100,000ppl all Ages, and split by age group
%loop through gridded TMean data to calculate daily mortality risk with
%constant and future population for future time periods.
   
avMortality = zeros(nRows,nGridCells);
avMortalitySE = zeros(nRows,nGridCells); %SE = Socioeconomic scenario i.e. population change
avMortPerIncrementUK = zeros(nClimScen*2,30); % results for CC only and CC+SE for plotting
header  = 1:30;
avMortPerIncrementUK(1,:) = header; % degrees above threshold in 1 degree intervals
h=1; %1:nGridCells

%% PAST
%Check if results already available - only run once then save.
if BiasCorrected == 0
    if exist ('ResultsPast.mat', 'file')
       load ('ResultsPast')
       disp ('loading Results: Past for plotting')
       
    elseif exist ('input_data_past.mat', 'file')  %calculate for each time period - based on input files
    disp ('Calculating mortality risk - Near Past')

        %Preallocate cell arrays
        ResultsPast = cell(adaptScen,1); %append all final gridded results for plotting and saving, for each age group, one array for each adaptation scenario
        for n = 1:adaptScen
            ResultsPast{n} = zeros(nGridCells, 17);
        end

        ResultsPastAgeGroup = cell(adaptScen,1); %append all final results aggregated per age group, for plotting and saving, one array for each adaptation scenario
        for n = 1:adaptScen
            ResultsPastAgeGroup{n} = zeros(ageGroupCnt, 3); 
        end

            % calculate mortality per gridcell
            for a = 1:adaptScen
                dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                for n = 1:nGridCells
                    dailyMort{n} = zeros(nRows, nCols);
                end

                totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                for n = 1:nGridCells
                    totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only        
                for n = 1:nGridCells
                    annMortality{n} = zeros(nRows, nCols);
                end
                
                for m = 1:ageGroupCnt %for all ages and per age group
                    i=1;
                    while (h <= nGridCells)
                        if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                            for r = 1:nRows
                                for c = 1:nCols
                                    if (data_past{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0) %plus 1 as impact calculated for each 1degreeC above defined threshold in baseline.
                                       dailyMort{h}(r,c) = ((((data_past{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1}(h,m+2); % Outputs the additional number of daily deaths %%increments 2 as population includes lat and lon coordinates. 
                                       incrementExceeded = round((data_past{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                       totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

        %% average annual mortality per increment %%for all ages only and no adaptation
                    if (m==1) && (a==1)      
                        for h = 1:nGridCells
                           for r = 1:nRows
                                avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                           end
                        end
                           dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                           matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                           meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells) %%change meanArray to sumArray
                           avMortPerIncrementUK(2,:) = mean(meanArray); % average across 12 RCMs
                    end

        %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                    for h = 1:nGridCells
                        annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).'; %re-shape array
                        %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                        for r = 1:nRows
                            avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                        end
                        %Absolute values
                        p_low = prctile(avMortality,(10),1);%10th percentile across years/GCMs
                        p_high = prctile(avMortality,(90),1);%90th percentile across years/GCMs
                        av = mean(avMortality,1);%average
                        sum_av = sum(av);
                        sum_p_low = sum(p_low);
                        sum_p_high = sum(p_high);
                    end
                        % Values per 100,000 ppl
                        av_per_population = sum_av/population_aggregate(m,1)*100000; %mean deaths per year/per 100,000 population - as a proportion of age group related population
                        p_low_per_population = sum_p_low/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 
                        p_high_per_population = sum_p_high/population_aggregate(m,1)*100000; %90th P deaths per year/per 1000 population 

                    %Save results for plotting - gridded
                    if m == 1 % all ages
                        ResultsPast{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                        ResultsPast{a}(1:h,1) = lonIDIndex(:,1); %longitude
                        ResultsPast{a}(1:h,2) = latIDIndex(:,1); %latitude
                    %Save results for plotting - UK aggregate
                        ResultsPastAgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];

                    elseif m==2
                        ResultsPast{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                        ResultsPastAgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                    elseif m==3
                        ResultsPast{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                        ResultsPastAgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                    elseif m==4
                        ResultsPast{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                        ResultsPastAgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                    elseif m==5
                        ResultsPast{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                        ResultsPastAgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                    end

                    %At this stage can also save daily data - don't use beyond calcs. above but other users
                    %may want daily data as output? For example...
        %             if (m == 1) && (a == 1) %all age groups, no adaptation
        %             S.dailyMortality = dailyMort;
        %             S.totalMortalityPerIndex = totMortPerIncrement;
        %             S.longitude = lonIDIndex;
        %             S.latitude = latIDIndex;
        %             save('dailyMort_AllAges_NoAdapt_CCOnly_Past.mat', '-struct', 'S','-v7.3')
        %             end

                    %Re-set counters
                    h=1; %1:nGridCells
                end
                    % Files to SAVE
                    save('ResultsPast.mat', 'ResultsPast', 'ResultsPastAgeGroup')
            end 
     end
   
    if exist ('Results_1_5C.mat', 'file')
       load ('Results_1_5C')
       load ('Results_1_5C_SE')
       disp ('loading Results: 1.5C for plotting')

    elseif exist ('input_data_1.5C.mat', 'file')
    disp ('Calculating mortality risk - 1.5DegreeC - CC Only')
        
       %set up cell array
        Results_1_5C = cell(adaptScen,1); %append all final results for plotting and saving, for each age group
        Results_1_5C_SE = cell(adaptScen,1);
        for n = 1:adaptScen
            Results_1_5C{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
            Results_1_5C_SE{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
        end

        Results_1_5C_AgeGroup = cell(adaptScen,1); %append all final results for plotting and saving, for each age group and new array for each adaptation scenario
        Results_1_5C_SE_AgeGroup = cell(adaptScen,1);
        for n = 1:adaptScen
            Results_1_5C_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
            Results_1_5C_SE_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
        end

        for l = 1:socioEcScen % Compute for CC-only and SE+CC
            if l == 1 %CC-only
                for a = 1:adaptScen
                    
                dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                for n = 1:nGridCells
                    dailyMort{n} = zeros(nRows, nCols);
                end

                totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                for n = 1:nGridCells
                    totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only        
                for n = 1:nGridCells
                    annMortality{n} = zeros(nRows, nCols);
                end

                    for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_1_5C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0)
                                           dailyMort{h}(r,c) = ((((data_1_5C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1}(h,m+2); % Outputs the additional number of daily deaths %%need to add population per gridcell and climate scenario in future
                                           incrementExceeded = round((data_1_5C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

                 %% average annual mortality per increment %%for all ages only
                        if (m==1) && (a==1)             
                            for h = 1:nGridCells
                               for r = 1:nRows
                                    avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                               end
                            end
                               dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                               matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                               meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                               avMortPerIncrementUK(3,:) = mean(meanArray);
                        end
                        
                    %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                        for h = 1:nGridCells
                            annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                            %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                            for r = 1:nRows
                                avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                            end
                            p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                            p_high = prctile(avMortality,(90),1);
                            av = mean(avMortality,1);
                            sum_av = sum(av);
                            sum_p_low = sum(p_low);
                            sum_p_high = sum(p_high);
                        end

                            av_per_population = sum_av/population_aggregate(m,1)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                            p_low_per_population = sum_p_low/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 
                            p_high_per_population = sum_p_high/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 

                        %Save results for plotting -
                        if m == 1
                            Results_1_5C{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                            Results_1_5C{a}(1:h,1) = lonIDIndex(:,1); %longitude
                            Results_1_5C{a}(1:h,2) = latIDIndex(:,1); %latitude

                            Results_1_5C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                        elseif m==2
                            Results_1_5C{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                            Results_1_5C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==3
                            Results_1_5C{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                            Results_1_5C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==4
                            Results_1_5C{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                            Results_1_5C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==5
                            Results_1_5C{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                            Results_1_5C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        end

                %Re-set counters
                        h=1; %1:nGridCells
                    end

                    % Files to SAVE
                    save('Results_1_5C.mat', 'Results_1_5C', 'Results_1_5C_AgeGroup')
                end
                
            elseif l == 2 %CC and Socio economic change (i.e. population change). Uses population data for 2020s
            disp ('Calculating mortality risk - 1.5DegreeC - CC+SE 2020s')

                for a = 1:adaptScen

                %Make sure arrays re-set to zero
                dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                for n = 1:nGridCells
                    dailyMort{n} = zeros(nRows, nCols);
                end

                totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                for n = 1:nGridCells
                    totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                end

                annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only        
                for n = 1:nGridCells
                    annMortality{n} = zeros(nRows, nCols);
                end

                    for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_1_5C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0)
                                           dailyMort{h}(r,c) = ((((data_1_5C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{2}(h,m+2); % Outputs the additional number of daily deaths
                                           incrementExceeded = round((data_1_5C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

                     %% average annual mortality per increment %%for all ages only
                        if (m==1) && (a==1)           
                            for h = 1:nGridCells
                               for r = 1:nRows
                                    avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                               end
                            end
                               dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                               matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                               meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                               avMortPerIncrementUK(4,:) = mean(meanArray);
                        end

                     %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                        for h = 1:nGridCells
                            annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                            %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                            for r = 1:nRows
                                avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                            end
                            p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                            p_high = prctile(avMortality,(90),1);
                            av = mean(avMortality,1);
                            sum_av = sum(av);
                            sum_p_low = sum(p_low);
                            sum_p_high = sum(p_high);
                        end

                            av_per_population = sum_av/population_aggregate(m,2)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                            p_low_per_population = sum_p_low/population_aggregate(m,2)*100000; %10th P deaths per year/per 1000 population 
                            p_high_per_population = sum_p_high/population_aggregate(m,2)*100000; %10th P deaths per year/per 1000 population 

                        %Save results for plotting -
                        if m == 1
                            Results_1_5C_SE{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                            Results_1_5C_SE{a}(1:h,1) = lonIDIndex(:,1); %longitude
                            Results_1_5C_SE{a}(1:h,2) = latIDIndex(:,1); %latitude

                            Results_1_5C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                        elseif m==2
                            Results_1_5C_SE{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                            Results_1_5C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==3
                            Results_1_5C_SE{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                            Results_1_5C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==4
                            Results_1_5C_SE{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                            Results_1_5C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==5
                            Results_1_5C_SE{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                            Results_1_5C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        end

                %Re-set counters
                        h=1; %1:nGridCells
                    end

                    % Files to SAVE
                    save('Results_1_5C_SE.mat', 'Results_1_5C_SE', 'Results_1_5C_SE_AgeGroup') 
                end

            end

        end
    end

    if exist ('Results_2C.mat' , 'file')
       load ('Results_2C')
       load ('Results_2C_SE')
       disp ('loading Results: 2.0C for plotting')

    elseif exist ('input_data_2.0C.mat', 'file')  %calculate for each time period - based on input files
    disp ('Calculating mortality risk - 2.0DegreeC - CC Only')

    % set up cell array
    Results_2C = cell(adaptScen,1); %append all final results for plotting and saving, for each age group
    Results_2C_SE = cell(adaptScen,1);
    for n = 1:adaptScen
        Results_2C{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
        Results_2C_SE{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
    end

    Results_2C_AgeGroup = cell(adaptScen,1); %append all final results for plotting and saving, for each age group and new array for each adaptation scenario
    Results_2C_SE_AgeGroup = cell(adaptScen,1);
    for n = 1:adaptScen
        Results_2C_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
        Results_2C_SE_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
    end

        for l = 1:socioEcScen % Compute for CC-only and SE+CC
            if l == 1 %CC-only
                for a = 1:adaptScen
                    
                    % Pre-allocate Cell Arrays to store intermediate results 1:adaptScen
                    dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                    for n = 1:nGridCells
                        dailyMort{n} = zeros(nRows, nCols);
                    end

                    totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                    for n = 1:nGridCells
                        totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                    end
                    
                    annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only      
                    for n = 1:nGridCells
                        annMortality{n} = zeros(nRows, nCols);
                    end
                    
                    avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                    for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                    end

                    for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_2_0C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0) 
                                           dailyMort{h}(r,c) = ((((data_2_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1}(h,m+2); % Outputs the additional number of daily deaths %%need to add population per gridcell and climate scenario in future
                                           incrementExceeded = round((data_2_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

        %% average annual mortality per increment %%for all ages only
                if (m==1) && (a==1)          
                    for h = 1:nGridCells
                       for r = 1:nRows
                            avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                       end
                    end
                       dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                       matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                       meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                       avMortPerIncrementUK(5,:) = mean(meanArray);
                end

            %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                for h = 1:nGridCells
                    annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                   % annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                    for r = 1:nRows
                        avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                    p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                    p_high = prctile(avMortality,(90),1);
                    av = mean(avMortality,1);
                    sum_av = sum(av);
                    sum_p_low = sum(p_low);
                    sum_p_high = sum(p_high);
                end

                    av_per_population = sum_av/population_aggregate(m,1)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                    p_low_per_population = sum_p_low/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 

                %Save results for plotting -
                if m == 1
                    Results_2C{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                    Results_2C{a}(1:h,1) = lonIDIndex(:,1); %longitude
                    Results_2C{a}(1:h,2) = latIDIndex(:,1); %latitude

                    Results_2C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                elseif m==2
                    Results_2C{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                    Results_2C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==3
                    Results_2C{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                    Results_2C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==4
                    Results_2C{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                    Results_2C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==5
                    Results_2C{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                    Results_2C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                end

                %Re-set counters
                        h=1; %1:nGridCells
                    end

                    % Files to SAVE
                    save('Results_2C.mat', 'Results_2C', 'Results_2C_AgeGroup')
                end

            elseif l == 2 %CC and Socio economic change (i.e. population change). Uses population data for 2020s
            disp ('Calculating mortality risk - 2.0DegreeC - CC+SE 2030s')

                for a = 1:adaptScen
                    % Pre-allocate Cell Arrays to store intermediate results 1:adaptScen
                    dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                    for n = 1:nGridCells
                        dailyMort{n} = zeros(nRows, nCols);
                    end
                    
                    annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only      
                    for n = 1:nGridCells
                        annMortality{n} = zeros(nRows, nCols);
                    end

                    totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                    for n = 1:nGridCells
                        totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                    end
                    
                    avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                    for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                    end

                    for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_2_0C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0)
                                           dailyMort{h}(r,c) = ((((data_2_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{3}(h,m+2); % Outputs the additional number of daily deaths %%need to add population per gridcell and climate scenario in future
                                           incrementExceeded = round((data_2_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

%% average annual mortality per increment %%for all ages only
                if (m==1) && (a==1)           
                    for h = 1:nGridCells
                       for r = 1:nRows
                            avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                       end
                    end
                       dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                       matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                       meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                       avMortPerIncrementUK(6,:) = mean(meanArray);
                end

            %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                for h = 1:nGridCells
                    annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                    %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                    for r = 1:nRows
                        avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                    p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                    p_high = prctile(avMortality,(90),1);
                    av = mean(avMortality,1);
                    sum_av = sum(av);
                    sum_p_low = sum(p_low);
                    sum_p_high = sum(p_high);
                end

                    av_per_population = sum_av/population_aggregate(m,3)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                    p_low_per_population = sum_p_low/population_aggregate(m,3)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,3)*100000; %10th P deaths per year/per 1000 population 

                %Save results for plotting -
                if m == 1
                    Results_2C_SE{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                    Results_2C_SE{a}(1:h,1) = lonIDIndex(:,1); %longitude
                    Results_2C_SE{a}(1:h,2) = latIDIndex(:,1); %latitude

                    Results_2C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                elseif m==2
                    Results_2C_SE{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                    Results_2C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==3
                    Results_2C_SE{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                    Results_2C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==4
                    Results_2C_SE{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                    Results_2C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==5
                    Results_2C_SE{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                    Results_2C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                end

                %Re-set counters
                        h=1; %1:nGridCells
                    end

                    % Files to SAVE
                    save('Results_2C_SE.mat', 'Results_2C_SE', 'Results_2C_SE_AgeGroup')
                end
            end
        end
    end

    if exist ('Results_3C.mat', 'file')
       load ('Results_3C')
       load ('Results_3C_SE')
       disp ('loading Results: 3.0C for plotting')

    elseif exist ('input_data_3.0C.mat', 'file')  %calculate for each time period - based on input files
    disp ('Calculating mortality risk - 3.0DegreeC - CC Only')

        %set up cell array
        Results_3C = cell(adaptScen,1); %append all final results for plotting and saving, for each age group
        Results_3C_SE = cell(adaptScen,1);
        for n = 1:adaptScen
            Results_3C{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
            Results_3C_SE{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
        end

        Results_3C_AgeGroup = cell(adaptScen,1); %append all final results for plotting and saving, for each age group and new array for each adaptation scenario
        Results_3C_SE_AgeGroup = cell(adaptScen,1);
        for n = 1:adaptScen
            Results_3C_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
            Results_3C_SE_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
        end

        for l = 1:socioEcScen % Compute for CC-only and SE+CC
            if l == 1 %CC-only
                for a = 1:adaptScen

                    dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                    for n = 1:nGridCells
                        dailyMort{n} = zeros(nRows, nCols);
                    end

                    totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                    for n = 1:nGridCells
                        totMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                    end

                    avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                    for n = 1:nGridCells
                        avMortPerIncrement{n} = zeros(nRows, 30); %default number columns assuming temp wont exceed this range...
                    end

                    annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only        
                    for n = 1:nGridCells
                        annMortality{n} = zeros(nRows, nCols);
                    end


                    for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_3_0C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0) 
                                           dailyMort{h}(r,c) = ((((data_3_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1}(h,m+2); % Outputs the additional number of daily deaths %%need to add population per gridcell and climate scenario in future
                                           incrementExceeded = round((data_3_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

%% average annual mortality per increment %%for all ages only
                if (m==1) && (a==1)            
                    for h = 1:nGridCells
                       for r = 1:nRows
                            avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                       end   
                    end
                       dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                       matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                       meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                       avMortPerIncrementUK(7,:) = mean(meanArray);
                end

            %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                for h = 1:nGridCells
                    annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                    %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                    for r = 1:nRows
                        avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                    p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                    p_high = prctile(avMortality,(90),1);
                    av = mean(avMortality,1);
                    sum_av = sum(av);
                    sum_p_low = sum(p_low);
                    sum_p_high = sum(p_high);
                end

                    av_per_population = sum_av/population_aggregate(m,1)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                    p_low_per_population = sum_p_low/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 

                %Save results for plotting -
                if m == 1
                    Results_3C{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                    Results_3C{a}(1:h,1) = lonIDIndex(:,1); %longitude
                    Results_3C{a}(1:h,2) = latIDIndex(:,1); %latitude

                    Results_3C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                elseif m==2
                    Results_3C{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                    Results_3C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==3
                    Results_3C{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                    Results_3C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==4
                    Results_3C{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                    Results_3C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==5
                    Results_3C{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                    Results_3C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                end

                %Re-set counters
                        h=1; %1:nGridCells
                    end

                    % Files to SAVE
                    save('Results_3C.mat', 'Results_3C', 'Results_3C_AgeGroup')
                end

            elseif l == 2 %CC and Socio economic change (i.e. population change). Uses population data for 2020s
            disp ('Calculating mortality risk - 3DegreeC - CC+SE 2030s')

                for a = 1:adaptScen
                    % Pre-allocate Cell Arrays to store intermediate results 1:adaptScen
                    dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                    for n = 1:nGridCells
                        dailyMort{n} = zeros(nRows, nCols);
                    end

                    totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                    for n = 1:nGridCells
                        totMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                    end

                    avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                    for n = 1:nGridCells
                        avMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                    end

                    annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only        
                    for n = 1:nGridCells
                        annMortality{n} = zeros(nRows, nCols);
                    end

                     for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_3_0C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0) 
                                           dailyMort{h}(r,c) = ((((data_3_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{4}(h,m+2); % Outputs the additional number of daily deaths %%need to add population per gridcell and climate scenario in future
                                           incrementExceeded = round((data_3_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

%% average annual mortality per increment %%for all ages only
                if (m==1) && (a==1)            
                    for h = 1:nGridCells
                       for r = 1:nRows
                            avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                       end   
                    end
                       dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                       matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                       meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                       avMortPerIncrementUK(8,:) = mean(meanArray);
                end

            %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                for h = 1:nGridCells
                    annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                    %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                    for r = 1:nRows
                        avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                    p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                    p_high = prctile(avMortality,(90),1);
                    av = mean(avMortality,1);
                    sum_av = sum(av);
                    sum_p_low = sum(p_low);
                    sum_p_high = sum(p_high);
                end

                    av_per_population = sum_av/population_aggregate(m,4)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                    p_low_per_population = sum_p_low/population_aggregate(m,4)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,4)*100000; %10th P deaths per year/per 1000 population 

                %Save results for plotting -
                if m == 1
                    Results_3C_SE{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                    Results_3C_SE{a}(1:h,1) = lonIDIndex(:,1); %longitude
                    Results_3C_SE{a}(1:h,2) = latIDIndex(:,1); %latitude

                    Results_3C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                elseif m==2
                    Results_3C_SE{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                    Results_3C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==3
                    Results_3C_SE{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                    Results_3C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==4
                    Results_3C_SE{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                    Results_3C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==5
                    Results_3C_SE{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                    Results_3C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                end


                %Re-set counters
                        h=1; %1:nGridCells
                    end

                    % Files to SAVE
                    save('Results_3C_SE.mat', 'Results_3C_SE', 'Results_3C_SE_AgeGroup')
                end
            end
        end
    end
    
save ('avMortPerIncrementUK.mat', 'avMortPerIncrementUK')

elseif BiasCorrected == 1
    if exist ('ResultsPast_Bias.mat', 'file')
       load ('ResultsPast')
       disp ('loading Results: Past (bias corrected) for plotting')
       
    elseif exist ('input_data_past_bias.mat', 'file')  %calculate for each time period - based on input files
    disp ('Calculating mortality risk - Near Past (Bias corrected)')
    
    %Preallocate cell arrays
    ResultsPast = cell(adaptScen,1); %append all final gridded results for plotting and saving, for each age group, one array for each adaptation scenario
    for n = 1:adaptScen
        ResultsPast{n} = zeros(nGridCells, 17);
    end

    ResultsPastAgeGroup = cell(adaptScen,1); %append all final results aggregated per age group, for plotting and saving, one array for each adaptation scenario
    for n = 1:adaptScen
        ResultsPastAgeGroup{n} = zeros(ageGroupCnt, 3); 
    end

        % calculate mortality per gridcell
        for a = 1:adaptScen

                dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                for n = 1:nGridCells
                    dailyMort{n} = zeros(nRows, nCols);
                end

                totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                for n = 1:nGridCells
                    totMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                end

                avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                end

                annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only        
                for n = 1:nGridCells
                    annMortality{n} = zeros(nRows, nCols);
                end

            for m = 1:ageGroupCnt %for all ages and per age group
                i=1;
                while (h <= nGridCells)
                    if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                        for r = 1:nRows
                            for c = 1:nCols
                                if (data_past{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0) %plus 1 as impact calculated for each 1degreeC above defined threshold in baseline.
                                   dailyMort{h}(r,c) = ((((data_past{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1}(h,m+2); % Outputs the additional number of daily deaths %%increments 2 as population includes lat and lon coordinates. 
                                   incrementExceeded = round((data_past{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                   totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

%% average annual mortality per increment %%for all ages only
                if (m==1) && (a==1)             
                    for h = 1:nGridCells
                       for r = 1:nRows
                            avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                       end
                    end
                       dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                       matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                       meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                       avMortPerIncrementUK(2,:) = mean(meanArray);
                end

    %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                for h = 1:nGridCells
                    annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                    %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                    for r = 1:nRows
                        avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                    %Absolute values
                    p_low = prctile(avMortality,(10),1);%10th percentile across years/GCMs
                    p_high = prctile(avMortality,(90),1);%90th percentile across years/GCMs
                    av = mean(avMortality,1);%average
                    sum_av = sum(av);
                    sum_p_low = sum(p_low);
                    sum_p_high = sum(p_high);
                end
                    % Values per 100,000 ppl
                    av_per_population = sum_av/population_aggregate(m,1)*100000; %mean deaths per year/per 100,000 population - as a proportion of age group related population
                    p_low_per_population = sum_p_low/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,1)*100000; %90th P deaths per year/per 1000 population 

                %Save results for plotting - gridded
                if m == 1 % all ages
                    ResultsPast{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                    ResultsPast{a}(1:h,1) = lonIDIndex(:,1); %longitude
                    ResultsPast{a}(1:h,2) = latIDIndex(:,1); %latitude
                %Save results for plotting - UK aggregate
                    ResultsPastAgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];

                elseif m==2
                    ResultsPast{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                    ResultsPastAgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==3
                    ResultsPast{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                    ResultsPastAgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==4
                    ResultsPast{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                    ResultsPastAgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==5
                    ResultsPast{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                    ResultsPastAgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                end

                %Re-set counters
                h=1; %1:nGridCells
            end
                % Files to SAVE
                save('ResultsPast_Bias.mat', 'ResultsPast', 'ResultsPastAgeGroup') 
        end
    end

    if exist ('Results_1_5C_Bias.mat', 'file')
       load ('Results_1_5C')
       load ('Results_1_5C_SE')
       disp ('loading Results: 1.5C for plotting (Bias corrrected)')

    elseif exist ('input_data_1.5C_bias.mat', 'file')
    disp ('Calculating mortality risk - 1.5DegreeC (Bias)- CC Only')
        
   %set up cell array
    Results_1_5C = cell(adaptScen,1); %append all final results for plotting and saving, for each age group
    Results_1_5C_SE = cell(adaptScen,1);
    for n = 1:adaptScen
        Results_1_5C{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
        Results_1_5C_SE{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
    end

    Results_1_5C_AgeGroup = cell(adaptScen,1); %append all final results for plotting and saving, for each age group and new array for each adaptation scenario
    Results_1_5C_SE_AgeGroup = cell(adaptScen,1);
    for n = 1:adaptScen
        Results_1_5C_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
        Results_1_5C_SE_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
    end

        for l = 1:socioEcScen % Compute for CC-only and SE+CC
            if l == 1 %CC-only
                for a = 1:adaptScen
                    % Pre-allocate Cell Arrays to store intermediate results 1:adaptScen
                dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                for n = 1:nGridCells
                    dailyMort{n} = zeros(nRows, nCols);
                end

                totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                for n = 1:nGridCells
                    totMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                end

                avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                end

                annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only        
                for n = 1:nGridCells
                    annMortality{n} = zeros(nRows, nCols);
                end

                    for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_1_5C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0)
                                           dailyMort{h}(r,c) = ((((data_1_5C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1}(h,m+2); % Outputs the additional number of daily deaths %%need to add population per gridcell and climate scenario in future
                                           incrementExceeded = round((data_1_5C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

    %% average annual mortality per increment %%for all ages only
                        if (m==1) && (a==1)             
                            for h = 1:nGridCells
                               for r = 1:nRows
                                    avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                               end
                            end
                               dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                               matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                               meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                               avMortPerIncrementUK(3,:) = mean(meanArray);
                        end

                    %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                        for h = 1:nGridCells
                            annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                            %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                            for r = 1:nRows
                                avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                            end
                            p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                            p_high = prctile(avMortality,(90),1);
                            av = mean(avMortality,1);
                            sum_av = sum(av);
                            sum_p_low = sum(p_low);
                            sum_p_high = sum(p_high);
                        end

                            av_per_population = sum_av/population_aggregate(m,1)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                            p_low_per_population = sum_p_low/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 
                            p_high_per_population = sum_p_high/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 

                        %Save results for plotting -
                        if m == 1
                            Results_1_5C{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                            Results_1_5C{a}(1:h,1) = lonIDIndex(:,1); %longitude
                            Results_1_5C{a}(1:h,2) = latIDIndex(:,1); %latitude

                            Results_1_5C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                        elseif m==2
                            Results_1_5C{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                            Results_1_5C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==3
                            Results_1_5C{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                            Results_1_5C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==4
                            Results_1_5C{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                            Results_1_5C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==5
                            Results_1_5C{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                            Results_1_5C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        end

                %Re-set counters
                        h=1; %1:nGridCells
                    end

                    % Files to SAVE
                    save('Results_1_5C_Bias.mat', 'Results_1_5C', 'Results_1_5C_AgeGroup') 
                end
                
            elseif l == 2 %CC and Socio economic change (i.e. population change). Uses population data for 2020s
            disp ('Calculating mortality risk - 1.5DegreeC (bias) - CC+SE 2020s')

                for a = 1:adaptScen
                % Pre-allocate Cell Arrays to store intermediate results 1:adaptScen
                dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                for n = 1:nGridCells
                    dailyMort{n} = zeros(nRows, nCols);
                end

                totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                for n = 1:nGridCells
                    totMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                end

                avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                end

                annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only        
                for n = 1:nGridCells
                    annMortality{n} = zeros(nRows, nCols);
                end

                for m = 1:ageGroupCnt %for all ages and per age group
                    i=1;
                    while (h <= nGridCells)
                        if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                            for r = 1:nRows
                                for c = 1:nCols
                                    if (data_1_5C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0)
                                       dailyMort{h}(r,c) = ((((data_1_5C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{2}(h,m+2); % Outputs the additional number of daily deaths
                                       incrementExceeded = round((data_1_5C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                       totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

%% average annual mortality per increment %%for all ages only
                    if (m==1) && (a==1)            
                        for h = 1:nGridCells
                           for r = 1:nRows
                                avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                           end
                        end
                           dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                           matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                           meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                           avMortPerIncrementUK(4,:) = mean(meanArray);
                    end

                 %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                    for h = 1:nGridCells
                        annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                        %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                        for r = 1:nRows
                            avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                        end
                        p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                        p_high = prctile(avMortality,(90),1);
                        av = mean(avMortality,1);
                        sum_av = sum(av);
                        sum_p_low = sum(p_low);
                        sum_p_high = sum(p_high);
                    end

                        av_per_population = sum_av/population_aggregate(m,2)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                        p_low_per_population = sum_p_low/population_aggregate(m,2)*100000; %10th P deaths per year/per 1000 population 
                        p_high_per_population = sum_p_high/population_aggregate(m,2)*100000; %10th P deaths per year/per 1000 population 

                    %Save results for plotting -
                    if m == 1
                        Results_1_5C_SE{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                        Results_1_5C_SE{a}(1:h,1) = lonIDIndex(:,1); %longitude
                        Results_1_5C_SE{a}(1:h,2) = latIDIndex(:,1); %latitude

                        Results_1_5C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                    elseif m==2
                        Results_1_5C_SE{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                        Results_1_5C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                    elseif m==3
                        Results_1_5C_SE{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                        Results_1_5C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                    elseif m==4
                        Results_1_5C_SE{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                        Results_1_5C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                    elseif m==5
                        Results_1_5C_SE{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                        Results_1_5C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                    end

            %Re-set counters
                    h=1; %1:nGridCells
                end

                    % Files to SAVE
                    save('Results_1_5C_SE_Bias.mat', 'Results_1_5C_SE', 'Results_1_5C_SE_AgeGroup') 
                end
            end
        end
    end
    
    if exist ('Results_2C_Bias.mat' , 'file')
       load ('Results_2C')
       load ('Results_2C_SE')
       disp ('loading Results: 2.0C for plotting')

    elseif exist ('input_data_2.0C.mat', 'file')  %calculate for each time period - based on input files
    disp ('Calculating mortality risk - 2.0DegreeC (bias) - CC Only')

    % set up cell array
    Results_2C = cell(adaptScen,1); %append all final results for plotting and saving, for each age group
    Results_2C_SE = cell(adaptScen,1);
    for n = 1:adaptScen
        Results_2C{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
        Results_2C_SE{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
    end

    Results_2C_AgeGroup = cell(adaptScen,1); %append all final results for plotting and saving, for each age group and new array for each adaptation scenario
    Results_2C_SE_AgeGroup = cell(adaptScen,1);
    for n = 1:adaptScen
        Results_2C_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
        Results_2C_SE_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
    end

        for l = 1:socioEcScen % Compute for CC-only and SE+CC
            if l == 1 %CC-only
                for a = 1:adaptScen
                    
                    % Pre-allocate Cell Arrays to store intermediate results 1:adaptScen
                    dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                    for n = 1:nGridCells
                        dailyMort{n} = zeros(nRows, nCols);
                    end

                    totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                    for n = 1:nGridCells
                        totMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                    end
                    
                    annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only      
                    for n = 1:nGridCells
                        annMortality{n} = zeros(nRows, nCols);
                    end
                    
                    avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                    for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                    end

                    for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_2_0C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0) 
                                           dailyMort{h}(r,c) = ((((data_2_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1}(h,m+2); % Outputs the additional number of daily deaths %%need to add population per gridcell and climate scenario in future
                                           incrementExceeded = round((data_2_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

    %% average annual mortality per increment %%for all ages only
                        if (m==1) && (a==1)             
                            for h = 1:nGridCells
                               for r = 1:nRows
                                    avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                               end
                            end
                               dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                               matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                               meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                               avMortPerIncrementUK(5,:) = mean(meanArray);
                        end

            %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                for h = 1:nGridCells
                    annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                   % annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                    for r = 1:nRows
                        avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                    p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                    p_high = prctile(avMortality,(90),1);
                    av = mean(avMortality,1);
                    sum_av = sum(av);
                    sum_p_low = sum(p_low);
                    sum_p_high = sum(p_high);
                end

                    av_per_population = sum_av/population_aggregate(m,1)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                    p_low_per_population = sum_p_low/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 

                %Save results for plotting -
                if m == 1
                    Results_2C{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                    Results_2C{a}(1:h,1) = lonIDIndex(:,1); %longitude
                    Results_2C{a}(1:h,2) = latIDIndex(:,1); %latitude

                    Results_2C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                elseif m==2
                    Results_2C{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                    Results_2C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==3
                    Results_2C{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                    Results_2C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==4
                    Results_2C{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                    Results_2C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==5
                    Results_2C{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                    Results_2C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                end

                %Re-set counters
                        h=1; %1:nGridCells
                    end

                    % Files to SAVE
                    save('Results_2C_Bias.mat', 'Results_2C', 'Results_2C_AgeGroup')
                end

            elseif l == 2 %CC and Socio economic change (i.e. population change). Uses population data for 2020s
            disp ('Calculating mortality risk - 2.0DegreeC (bias) - CC+SE 2030s')

                for a = 1:adaptScen
                    % Pre-allocate Cell Arrays to store intermediate results 1:adaptScen
                    dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                    for n = 1:nGridCells
                        dailyMort{n} = zeros(nRows, nCols);
                    end
                    
                    annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only      
                    for n = 1:nGridCells
                        annMortality{n} = zeros(nRows, nCols);
                    end

                    totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                    for n = 1:nGridCells
                        totMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                    end
                    
                    avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                    for n = 1:nGridCells
                    avMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                    end

                    for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_2_0C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0)
                                           dailyMort{h}(r,c) = ((((data_2_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{3}(h,m+2); % Outputs the additional number of daily deaths %%need to add population per gridcell and climate scenario in future
                                           incrementExceeded = round((data_2_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

    %% average annual mortality per increment %%for all ages only
                        if (m==1) && (a==1)            
                            for h = 1:nGridCells
                               for r = 1:nRows
                                    avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                               end
                            end
                               dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                               matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                               meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                               avMortPerIncrementUK(6,:) = mean(meanArray);
                        end

            %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                        for h = 1:nGridCells
                            annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                            %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                            for r = 1:nRows
                                avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                            end
                            p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                            p_high = prctile(avMortality,(90),1);
                            av = mean(avMortality,1);
                            sum_av = sum(av);
                            sum_p_low = sum(p_low);
                            sum_p_high = sum(p_high);
                        end

                            av_per_population = sum_av/population_aggregate(m,3)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                            p_low_per_population = sum_p_low/population_aggregate(m,3)*100000; %10th P deaths per year/per 1000 population 
                            p_high_per_population = sum_p_high/population_aggregate(m,3)*100000; %10th P deaths per year/per 1000 population 

                        %Save results for plotting -
                        if m == 1
                            Results_2C_SE{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                            Results_2C_SE{a}(1:h,1) = lonIDIndex(:,1); %longitude
                            Results_2C_SE{a}(1:h,2) = latIDIndex(:,1); %latitude

                            Results_2C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                        elseif m==2
                            Results_2C_SE{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                            Results_2C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==3
                            Results_2C_SE{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                            Results_2C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==4
                            Results_2C_SE{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                            Results_2C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        elseif m==5
                            Results_2C_SE{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                            Results_2C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                        end

                        %Re-set counters
                                h=1; %1:nGridCells
                    end

                    % Files to SAVE
                    save('Results_2C_SE_Bias.mat', 'Results_2C_SE', 'Results_2C_SE_AgeGroup')
                end
            end
        end
    end

    if exist ('Results_3C_Bias.mat', 'file')
       load ('Results_3C')
       load ('Results_3C_SE')
       disp ('loading Results: 3.0C for plotting')

    elseif exist ('input_data_3.0C.mat', 'file')  %calculate for each time period - based on input files
    disp ('Calculating mortality risk - 3.0DegreeC (bias) - CC Only')

        %set up cell array
        Results_3C = cell(adaptScen,1); %append all final results for plotting and saving, for each age group
        Results_3C_SE = cell(adaptScen,1);
        for n = 1:adaptScen
            Results_3C{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
            Results_3C_SE{n} = zeros(nGridCells, 17); %set automatically in future based on gridcell rows
        end

        Results_3C_AgeGroup = cell(adaptScen,1); %append all final results for plotting and saving, for each age group and new array for each adaptation scenario
        Results_3C_SE_AgeGroup = cell(adaptScen,1);
        for n = 1:adaptScen
            Results_3C_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
            Results_3C_SE_AgeGroup{n} = zeros(ageGroupCnt, 3); %set automatically in future based on gridcell rows
        end

        for l = 1:socioEcScen % Compute for CC-only and SE+CC
            if l == 1 %CC-only
                for a = 1:adaptScen

                    dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                    for n = 1:nGridCells
                        dailyMort{n} = zeros(nRows, nCols);
                    end

                    totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                    for n = 1:nGridCells
                        totMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                    end

                    avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                    for n = 1:nGridCells
                        avMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                    end

                    annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only        
                    for n = 1:nGridCells
                        annMortality{n} = zeros(nRows, nCols);
                    end

                    for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_3_0C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0) 
                                           dailyMort{h}(r,c) = ((((data_3_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{1}(h,m+2); % Outputs the additional number of daily deaths %%need to add population per gridcell and climate scenario in future
                                           incrementExceeded = round((data_3_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

    %% average annual mortality per increment %%for all ages only
                        if (m==1) && (a==1)             
                            for h = 1:nGridCells
                               for r = 1:nRows
                                    avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                               end
                            end
                               dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                               matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                               meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                               avMortPerIncrementUK(7,:) = mean(meanArray);
                        end

            %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                for h = 1:nGridCells
                    annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                    %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                    for r = 1:nRows
                        avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                    p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                    p_high = prctile(avMortality,(90),1);
                    av = mean(avMortality,1);
                    sum_av = sum(av);
                    sum_p_low = sum(p_low);
                    sum_p_high = sum(p_high);
                end

                    av_per_population = sum_av/population_aggregate(m,1)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                    p_low_per_population = sum_p_low/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,1)*100000; %10th P deaths per year/per 1000 population 

                %Save results for plotting -
                if m == 1
                    Results_3C{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                    Results_3C{a}(1:h,1) = lonIDIndex(:,1); %longitude
                    Results_3C{a}(1:h,2) = latIDIndex(:,1); %latitude

                    Results_3C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                elseif m==2
                    Results_3C{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                    Results_3C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==3
                    Results_3C{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                    Results_3C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==4
                    Results_3C{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                    Results_3C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==5
                    Results_3C{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                    Results_3C_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                end

                %Re-set counters
                        h=1; %1:nGridCells
                        %i=1; %1:nGridCells
                        %k=1;
                    end

                    % Files to SAVE
                    save('Results_3C_Bias.mat', 'Results_3C', 'Results_3C_AgeGroup')
                end

            elseif l == 2 %CC and Socio economic change (i.e. population change). Uses population data for 2020s
            disp ('Calculating mortality risk - 3DegreeC (bias) - CC+SE 2030s')

                for a = 1:adaptScen
                    % Pre-allocate Cell Arrays to store intermediate results 1:adaptScen
                    dailyMort = cell(nGridCells,1);   %store intermediate value of daily additional deaths CC-only          
                    for n = 1:nGridCells
                        dailyMort{n} = zeros(nRows, nCols);
                    end

                    totMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (summed across all years/GCMs)CC-only              
                    for n = 1:nGridCells
                        totMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                    end

                    avMortPerIncrement = cell(nGridCells,1);   %store intermediate value of additional deaths split by each increment exceeded (average annual) CC only    
                    for n = 1:nGridCells
                        avMortPerIncrement{n} = zeros(nRows+1, 30); %default number columns assuming temp wont exceed this range...
                    end

                    annMortality = cell(nGridCells,1);   %store intermediate value of Mortality risk * population (average annual)CC only        
                    for n = 1:nGridCells
                        annMortality{n} = zeros(nRows, nCols);
                    end

                    for m = 1:ageGroupCnt %for all ages and per age group
                        i=1;
                        while (h <= nGridCells)
                            if gridID(i,1) == lonIDIndex(h,1) && gridID(i,2) == latIDIndex(h,1) % could just set dataup so matches but this looks through cellID to match to data_past as some grid cells may be missing, or may be in differet order in future files etc.
                                for r = 1:nRows
                                    for c = 1:nCols
                                        if (data_3_0C{i,2}(r,c) >= Mortality_Threshold_RR(h,4) + 1 + adaptIncrement(a)) && (Mortality_Threshold_RR(h,4) > 0)
                                           dailyMort{h}(r,c) = ((((data_3_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4))*Mortality_Threshold_RR(h,m+4))/100)*dailyDeathRate(h,m+2))*population{4}(h,m+2); % Outputs the additional number of daily deaths %%need to add population per gridcell and climate scenario in future
                                           incrementExceeded = round((data_3_0C{i,2}(r,c)- Mortality_Threshold_RR(h,4)));
                                           totMortPerIncrement{h}(r, incrementExceeded) = totMortPerIncrement{h}(r, incrementExceeded) + dailyMort{h}(r,c); %cumulative mortality per each 1 degree increment above the threshold
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

    %% average annual mortality per increment %%for all ages only
                        if (m==1) && (a==1)           
                            for h = 1:nGridCells
                               for r = 1:nRows
                                    avMortPerIncrement{h}(r,:) = totMortPerIncrement{h}(r,:)/years;
                               end
                            end
                               dim = ndims(avMortPerIncrement{1});          % Get the number of dimensions for arrays
                               matrix = cat(dim+1,avMortPerIncrement{:});        % Convert to a (dim+1)-dimensional matrix
                               meanArray = sum(matrix,dim+1);  % Get the mean across 3d array (summing mortality across UK grid cells)
                               avMortPerIncrementUK(8,:) = mean(meanArray);
                        end

            %% Calculate Annual risk (across years and across 12 GCMs) and percentile range
                for h = 1:nGridCells
                    annMortality{h} = squeeze(nansum(reshape(dailyMort{h}.',daysPerYear,size(dailyMort{h},2)./daysPerYear,[]))).';
                    %annMortalityDistribution{m}(:,:) = annMortalityDistribution{m}(:,:) + annMortality{h}(:,:); % used to look at distribution of all years
                    for r = 1:nRows
                        avMortality(r,h) = mean(annMortality{h}(r,:),[1,2]); %mean across years for each GCM
                    end
                    p_low = prctile(avMortality,(10),1);%10th and 90th percentile across years/GCMs
                    p_high = prctile(avMortality,(90),1);
                    av = mean(avMortality,1);
                    sum_av = sum(av);
                    sum_p_low = sum(p_low);
                    sum_p_high = sum(p_high);
                end

                    av_per_population = sum_av/population_aggregate(m,4)*100000; %mean deaths per year/per 100,000 population - as a proportion of total population
                    p_low_per_population = sum_p_low/population_aggregate(m,4)*100000; %10th P deaths per year/per 1000 population 
                    p_high_per_population = sum_p_high/population_aggregate(m,4)*100000; %10th P deaths per year/per 1000 population 

                %Save results for plotting -
                if m == 1
                    Results_3C_SE{a}(1:h,3:5) = [av.', p_low.', p_high.']; %3 variables[nGridCells*3] * 5 Age_group categories (m)
                    Results_3C_SE{a}(1:h,1) = lonIDIndex(:,1); %longitude
                    Results_3C_SE{a}(1:h,2) = latIDIndex(:,1); %latitude

                    Results_3C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.']; % UK aggregate per 100,000

                elseif m==2
                    Results_3C_SE{a}(1:h,6:8) = [av.', p_low.', p_high.'];
                    Results_3C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==3
                    Results_3C_SE{a}(1:h,9:11) = [av.', p_low.', p_high.'];
                    Results_3C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==4
                    Results_3C_SE{a}(1:h,12:14) = [av.', p_low.', p_high.'];
                    Results_3C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                elseif m==5
                    Results_3C_SE{a}(1:h,15:17) = [av.', p_low.', p_high.'];
                    Results_3C_SE_AgeGroup{a}(m, 1:3) = [av_per_population.', p_low_per_population.', p_high_per_population.'];
                end

                %Re-set counters
                        h=1; %1:nGridCells
                    end

                    % Files to SAVE
                    save('Results_3C_SE_Bias.mat', 'Results_3C_SE', 'Results_3C_SE_AgeGroup')
                end
            end
        end
    end
save ('avMortPerIncrementUK_bias.mat', 'avMortPerIncrementUK')
end

%% Amalgamate results for plotting - UK data based on sum of gridded data (All ages and by age group)Average, 10th and 90th P summed.
%Absolute values also split by age group - differs to mortality per 100,000
%population per age group.

AvAnnualDeathsUK = zeros(3,12);
AvAnnualDeathsUK_SE = zeros(3,12);
AvAnnualDeathsAdapt = zeros(3,4);
AvAnnualDeathsAdapt_SE = zeros(3,4);

for a = 1  %% update for adaptation scenarios in future....
        % All ages
    AvAnnualDeathsUK(1,1:3) = [sum(ResultsPast{a}(:,3)), sum(ResultsPast{a}(:,4)), sum(ResultsPast{a}(:,5))]; % Near Past
    AvAnnualDeathsUK(1,4:6) = [sum(Results_1_5C{a}(:,3)), sum(Results_1_5C{a}(:,4)), sum(Results_1_5C{a}(:,5))]; % 1.5 degree
    AvAnnualDeathsUK(1,7:9) = [sum(Results_2C{a}(:,3)), sum(Results_2C{a}(:,4)), sum(Results_2C{a}(:,5))]; % 2 degree C
    AvAnnualDeathsUK(1,10:12) = [sum(Results_3C{a}(:,3)), sum(Results_3C{a}(:,4)), sum(Results_3C{a}(:,5))]; % 3 degree C
    AvAnnualDeathsUK_SE(1,1:3) = [sum(ResultsPast{a}(:,3)), sum(ResultsPast{a}(:,4)), sum(ResultsPast{a}(:,5))];% Near Past (population same in both scenarios)
    AvAnnualDeathsUK_SE(1,4:6) = [sum(Results_1_5C_SE{a}(:,3)), sum(Results_1_5C_SE{a}(:,4)), sum(Results_1_5C_SE{a}(:,5))]; % 1.5 degree CC+SE
    AvAnnualDeathsUK_SE(1,7:9) = [sum(Results_2C_SE{a}(:,3)), sum(Results_2C_SE{a}(:,4)), sum(Results_2C_SE{a}(:,5))]; %2 degree C CC+SE
    AvAnnualDeathsUK_SE(1,10:12) = [sum(Results_3C_SE{a}(:,3)), sum(Results_3C_SE{a}(:,4)), sum(Results_3C_SE{a}(:,5))]; % 3 degree C CC+SE

    %To calculate absolute deaths per age group as opposed to per 100,000
    %ppl
    AvAnnualDeathsUK_age = zeros(4,4); %set up arrays
    AvAnnualDeathsUK_SE_age = zeros(4,4);
    AvAnnualDeathsUK_age_10th = zeros(4,4);
    AvAnnualDeathsUK_SE_age_10th = zeros(4,4);
    AvAnnualDeathsUK_age_90th = zeros(4,4);
    AvAnnualDeathsUK_SE_age_90th = zeros(4,4);

    AvAnnualDeathsUK_age(1,1:4) = [sum(ResultsPast{a}(:,6)), sum(ResultsPast{a}(:,9)), sum(ResultsPast{a}(:,12)), sum(ResultsPast{a}(:,15))]; % Near Past
    AvAnnualDeathsUK_age(2,1:4) =[sum(Results_1_5C{a}(:,6)), sum(Results_1_5C{a}(:,9)), sum(Results_1_5C{a}(:,12)), sum(Results_1_5C{a}(:,15))]; % 1.5 degree
    AvAnnualDeathsUK_age(3,1:4) = [sum(Results_2C{a}(:,6)), sum(Results_2C{a}(:,9)), sum(Results_2C{a}(:,12)), sum(Results_2C{a}(:,15))];  % 2 degree C
    AvAnnualDeathsUK_age(4,1:4) = [sum(Results_3C{a}(:,6)), sum(Results_3C{a}(:,9)), sum(Results_3C{a}(:,12)), sum(Results_3C{a}(:,15))]; % 3 degree C
    AvAnnualDeathsUK_SE_age(1,1:4) = AvAnnualDeathsUK_age(1,1:4);% Near Past (population same in both scenarios)
    AvAnnualDeathsUK_SE_age(2,1:4) = [sum(Results_1_5C_SE{a}(:,6)), sum(Results_1_5C_SE{a}(:,9)), sum(Results_1_5C_SE{a}(:,12)),sum(Results_1_5C_SE{a}(:,15))]; % 1.5 degree CC+SE
    AvAnnualDeathsUK_SE_age(3,1:4) = [sum(Results_2C_SE{a}(:,6)), sum(Results_2C_SE{a}(:,9)), sum(Results_2C_SE{a}(:,12)),sum(Results_2C_SE{a}(:,15))]; %2 degree C CC+SE
    AvAnnualDeathsUK_SE_age(4,1:4) = [sum(Results_3C_SE{a}(:,6)), sum(Results_3C_SE{a}(:,9)), sum(Results_3C_SE{a}(:,12)),sum(Results_3C_SE{a}(:,15))]; % 3 degree C CC+SE

    AvAnnualDeathsUK_age_10th(1,1:4) = [sum(ResultsPast{a}(:,7)), sum(ResultsPast{a}(:,10)), sum(ResultsPast{a}(:,13)), sum(ResultsPast{a}(:,16))]; % Near Past
    AvAnnualDeathsUK_age_10th(2,1:4) =[sum(Results_1_5C{a}(:,7)), sum(Results_1_5C{a}(:,10)), sum(Results_1_5C{a}(:,13)), sum(Results_1_5C{a}(:,16))]; % 1.5 degree
    AvAnnualDeathsUK_age_10th(3,1:4) = [sum(Results_2C{a}(:,7)), sum(Results_2C{a}(:,10)), sum(Results_2C{a}(:,13)), sum(Results_2C{a}(:,16))];  % 2 degree C
    AvAnnualDeathsUK_age_10th(4,1:4) = [sum(Results_3C{a}(:,7)), sum(Results_3C{a}(:,10)), sum(Results_3C{a}(:,13)), sum(Results_3C{a}(:,16))]; % 3 degree Csum(Results_2C{a}(:,15));
    AvAnnualDeathsUK_SE_age_10th(1,1:4) = AvAnnualDeathsUK_age_10th(1,1:4);% Near Past (population same in both scenarios)
    AvAnnualDeathsUK_SE_age_10th(2,1:4) = [sum(Results_1_5C_SE{a}(:,7)), sum(Results_1_5C_SE{a}(:,10)), sum(Results_1_5C_SE{a}(:,13)),sum(Results_1_5C_SE{a}(:,16))]; % 1.5 degree CC+SE
    AvAnnualDeathsUK_SE_age_10th(3,1:4) = [sum(Results_2C_SE{a}(:,7)), sum(Results_2C_SE{a}(:,10)), sum(Results_2C_SE{a}(:,13)),sum(Results_2C_SE{a}(:,16))]; %2 degree C CC+SE
    AvAnnualDeathsUK_SE_age_10th(4,1:4) = [sum(Results_3C_SE{a}(:,7)), sum(Results_3C_SE{a}(:,10)), sum(Results_3C_SE{a}(:,13)),sum(Results_3C_SE{a}(:,16))]; % 3 degree C CC+SE

    AvAnnualDeathsUK_age_90th(1,1:4) = [sum(ResultsPast{a}(:,8)), sum(ResultsPast{a}(:,11)), sum(ResultsPast{a}(:,14)), sum(ResultsPast{a}(:,17))]; % Near Past
    AvAnnualDeathsUK_age_90th(2,1:4) =[sum(Results_1_5C{a}(:,8)), sum(Results_1_5C{a}(:,11)), sum(Results_1_5C{a}(:,14)), sum(Results_1_5C{a}(:,17))]; % 1.5 degree
    AvAnnualDeathsUK_age_90th(3,1:4) = [sum(Results_2C{a}(:,8)), sum(Results_2C{a}(:,11)), sum(Results_2C{a}(:,14)), sum(Results_2C{a}(:,17))];  % 2 degree C
    AvAnnualDeathsUK_age_90th(4,1:4) = [sum(Results_3C{a}(:,8)), sum(Results_3C{a}(:,11)), sum(Results_3C{a}(:,14)), sum(Results_3C{a}(:,17))]; % 3 degree Csum(Results_2C{a}(:,15));
    AvAnnualDeathsUK_SE_age_90th(1,1:4) = AvAnnualDeathsUK_age_90th(1,1:4);% Near Past (population same in both scenarios)
    AvAnnualDeathsUK_SE_age_90th(2,1:4) = [sum(Results_1_5C_SE{a}(:,8)), sum(Results_1_5C_SE{a}(:,11)), sum(Results_1_5C_SE{a}(:,14)),sum(Results_1_5C_SE{a}(:,17))]; % 1.5 degree CC+SE
    AvAnnualDeathsUK_SE_age_90th(3,1:4) = [sum(Results_2C_SE{a}(:,8)), sum(Results_2C_SE{a}(:,11)), sum(Results_2C_SE{a}(:,14)),sum(Results_2C_SE{a}(:,17))]; %2 degree C CC+SE
    AvAnnualDeathsUK_SE_age_90th(4,1:4) = [sum(Results_3C_SE{a}(:,8)), sum(Results_3C_SE{a}(:,11)), sum(Results_3C_SE{a}(:,14)),sum(Results_3C_SE{a}(:,17))]; % 3 degree C CC+SE
end
    %difference between no adaptation and adaptation strategies for stacked
    %plot
    AvAnnualDeathsAdapt(1,1:4) = [sum(ResultsPast{1}(:,3))-sum(ResultsPast{2}(:,3)), sum(Results_1_5C{1}(:,3))-sum(Results_1_5C{2}(:,3)), sum(Results_2C{1}(:,3))-sum(Results_2C{2}(:,3)), sum(Results_3C{1}(:,3)) - sum(Results_3C{2}(:,3))]; 
    AvAnnualDeathsAdapt(2,1:4) = [sum(ResultsPast{2}(:,3))- sum(ResultsPast{3}(:,3)), sum(Results_1_5C{2}(:,3))-sum(Results_1_5C{3}(:,3)), sum(Results_2C{2}(:,3))-sum(Results_2C{3}(:,3)), sum(Results_3C{2}(:,3))-sum(Results_3C{3}(:,3))];
    AvAnnualDeathsAdapt(3,1:4) = [sum(ResultsPast{3}(:,3)), sum(Results_1_5C{3}(:,3)), sum(Results_2C{3}(:,3)), sum(Results_3C{3}(:,3))]; 

    AvAnnualDeathsAdapt_SE(1,1:4) = [sum(ResultsPast{1}(:,3))-sum(ResultsPast{2}(:,3)), sum(Results_1_5C_SE{1}(:,3))-sum(Results_1_5C_SE{2}(:,3)), sum(Results_2C_SE{1}(:,3))-sum(Results_2C_SE{2}(:,3)), sum(Results_3C_SE{1}(:,3))-sum(Results_3C_SE{2}(:,3))];
    AvAnnualDeathsAdapt_SE(2,1:4) = [sum(ResultsPast{2}(:,3))-sum(ResultsPast{3}(:,3)), sum(Results_1_5C_SE{2}(:,3))-sum(Results_1_5C_SE{3}(:,3)), sum(Results_2C_SE{2}(:,3))-sum(Results_2C_SE{3}(:,3)), sum(Results_3C_SE{2}(:,3))-sum(Results_3C_SE{3}(:,3))];
    AvAnnualDeathsAdapt_SE(3,1:4) = [sum(ResultsPast{3}(:,3)), sum(Results_1_5C_SE{3}(:,3)), sum(Results_2C_SE{3}(:,3)), sum(Results_3C_SE{3}(:,3))];
    
%% PLOTS
% Create some initial output figures

%1. Average annual deaths for each scenario (All ages) CC-only and CC+SE

figure;
M = max(AvAnnualDeathsUK_SE(1,:)*1.1); %set y axis max limit across sub-plots
x=1:nClimScen;
subplot(1,2,1);

% Create bar and set individual colours
b = bar(AvAnnualDeathsUK(1,[1,4,7,10]), 'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.843137264251709 0];
b.CData(2,:) = [0.87058824300766 0.490196079015732 0];
b.CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
b.CData(4,:) = [0.600000023841858 0.200000002980232 0];

hold on
p_errhigh = AvAnnualDeathsUK(1,[3,6,9,12])- AvAnnualDeathsUK(1,[1,4,7,10]); %90th p - average
p_errlow = AvAnnualDeathsUK(1,[1,4,7,10])-AvAnnualDeathsUK(1,[2,5,8,11]); % average - 10th p)

er = errorbar(x, AvAnnualDeathsUK(1,[1,4,7,10]), p_errlow, p_errhigh, 'LineWidth',1);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

% Create ylabel
ylabel('Average Annual Heat Related Deaths');

% Create xlabel
xlabel('Climate Scenario');

% Create title
title('Climate Change Only');

% Set the remaining axes properties
set (gca, 'ylim', [0 M])
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Past','1.5','2.0','3.0'},'YGrid',...
    'on');

subplot(1,2,2);
% Create bar
b= bar(AvAnnualDeathsUK_SE(1,[1,4,7,10]), 'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.843137264251709 0];
b.CData(2,:) = [0.87058824300766 0.490196079015732 0];
b.CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
b.CData(4,:) = [0.600000023841858 0.200000002980232 0];

hold on
p_errhighSE = AvAnnualDeathsUK_SE(1,[3,6,9,12])- AvAnnualDeathsUK_SE(1,[1,4,7,10]); %90th p - average
p_errlowSE = AvAnnualDeathsUK_SE(1,[1,4,7,10])-AvAnnualDeathsUK_SE(1,[2,5,8,11]); % average - 10th p)

er = errorbar(x, AvAnnualDeathsUK_SE(1,[1,4,7,10]), p_errlowSE, p_errhighSE, 'LineWidth',1);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

% Create ylabel
%ylabel('Average Annual Heat Related Deaths');

% Create xlabel
xlabel('Climate Scenario');

% Create title
title('Climate and Population Change');

% Set the remaining axes properties
set (gca, 'ylim', [0 M])
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Past','1.5','2.0','3.0'},'YGrid',...
    'on');

%% Clustered bars per age group, climate scenario (cc and cc+SE same as weighted to population)
% Per 100,000 population

AvAgeResults = zeros(nClimScen, ageGroupCnt-1); %ignore all ages.
AvAgeResults_10percentile = zeros(nClimScen, ageGroupCnt-1);
AvAgeResults_90percentile = zeros(nClimScen, ageGroupCnt-1);
% AvAgeResultsSE = zeros(nClimScen, ageGroupCnt-1); %ignore all ages.
% AvAgeResults_10percetileSE = zeros(nClimScen, ageGroupCnt-1);
% AvAgeResults_90percetileSE = zeros(nClimScen, ageGroupCnt-1);

AvAgeResults(1,:) = ResultsPastAgeGroup{1}(2:5,1);
AvAgeResults(2,:) = Results_1_5C_AgeGroup{1}(2:5,1);
AvAgeResults(3,:) = Results_2C_AgeGroup{1}(2:5,1);
AvAgeResults(4,:) = Results_3C_AgeGroup{1}(2:5,1);

AvAgeResults_10percentile(1,:) = ResultsPastAgeGroup{1}(2:5,2);
AvAgeResults_10percentile(2,:) = Results_1_5C_AgeGroup{1}(2:5,2); 
AvAgeResults_10percentile(3,:) = Results_2C_AgeGroup{1}(2:5,2); 
AvAgeResults_10percentile(4,:) = Results_3C_AgeGroup{1}(2:5,2);

AvAgeResults_90percentile(1,:) = ResultsPastAgeGroup{1}(2:5,3);
AvAgeResults_90percentile(2,:) = Results_1_5C_AgeGroup{1}(2:5,3); 
AvAgeResults_90percentile(3,:) = Results_2C_AgeGroup{1}(2:5,3); 
AvAgeResults_90percentile(4,:) = Results_3C_AgeGroup{1}(2:5,3); 
% 
% AvAgeResultsSE(1,:) = ResultsPastAgeGroup{1}(2:5,1);
% AvAgeResultsSE(2,:) = Results_1_5C_SE_AgeGroup{1}(2:5,1);
% AvAgeResultsSE(3,:) = Results_2C_SE_AgeGroup{1}(2:5,1);
% AvAgeResultsSE(4,:) = Results_3C_SE_AgeGroup{1}(2:5,1);
% 
% AvAgeResults_10percentileSE(1,:) = ResultsPastAgeGroup{1}(2:5,2); 
% AvAgeResults_10percentileSE(2,:) = Results_1_5C_SE_AgeGroup{1}(2:5,2); 
% AvAgeResults_10percentileSE(3,:) = Results_2C_SE_AgeGroup{1}(2:5,2); 
% AvAgeResults_10percentileSE(4,:) = Results_3C_SE_AgeGroup{1}(2:5,2);
% 
% AvAgeResults_90percentileSE(1,:) = ResultsPastAgeGroup{1}(2:5,3); 
% AvAgeResults_90percentileSE(2,:) = Results_1_5C_SE_AgeGroup{1}(2:5,3); 
% AvAgeResults_90percentileSE(3,:) = Results_2C_SE_AgeGroup{1}(2:5,3); 
% AvAgeResults_90percentileSE(4,:) = Results_3C_SE_AgeGroup{1}(2:5,3);

AvAgeResults = AvAgeResults.'; %Transpose for clustered bar
% AvAgeResultsSE = AvAgeResultsSE.'; %Transpose for clustered bar
AvAgeResults_10percentile = AvAgeResults_10percentile.';
AvAgeResults_90percentile = AvAgeResults_90percentile.';
% AvAgeResults_10percentileSE = AvAgeResults_10percentileSE.';
% AvAgeResults_90percentileSE = AvAgeResults_90percentileSE.';

M = max(AvAgeResults_90percentile(:,4)*1.1); %set y axis max limit across sub-plots
figure
%subplot(1,2,1);
b = bar(AvAgeResults, 'FaceColor', 'flat');
b(1).CData(1,:) = [1 0.843137264251709 0];
b(2).CData(1,:) = [0.87058824300766 0.490196079015732 0];
b(3).CData(1,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
b(4).CData(1,:) = [0.600000023841858 0.200000002980232 0];
b(1).CData(2,:) = [1 0.843137264251709 0];
b(2).CData(2,:) = [0.87058824300766 0.490196079015732 0];
b(3).CData(2,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
b(4).CData(2,:) = [0.600000023841858 0.200000002980232 0];
b(1).CData(3,:) = [1 0.843137264251709 0];
b(2).CData(3,:) = [0.87058824300766 0.490196079015732 0];
b(3).CData(3,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
b(4).CData(3,:) = [0.600000023841858 0.200000002980232 0];
b(1).CData(4,:) = [1 0.843137264251709 0];
b(2).CData(4,:) = [0.87058824300766 0.490196079015732 0];
b(3).CData(4,:) = [0.850980401039124 0.325490206480026 0.0980392172932625];
b(4).CData(4,:) = [0.600000023841858 0.200000002980232 0];

hold on
p_errhigh = AvAgeResults_90percentile - AvAgeResults; %90th p - average
p_errlow = AvAgeResults - AvAgeResults_10percentile; % average - 10th p)

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(AvAgeResults);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

er = errorbar(x', AvAgeResults, p_errlow, p_errhigh, 'k', 'LineWidth',1,'linestyle','none');                                

% Create ylabel
ylabel('Average Annual Heat Related Deaths per 100,000 population');


% Create xlabel
xlabel('Age Group');
%create title
%title ('Climate Change Only')


set (gca, 'ylim', [0 M])
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'0-64', '65-74', '75-84','85+'},'YGrid',...
    'on');

% subplot(1,2,2);
% b= bar(AvAgeResultsSE.'); %create bar
% hold on
% 
% p_errhigh = AvAgeResults_90percentileSE - AvAgeResultsSE; %90th p - average
% p_errlow = AvAgeResultsSE - AvAgeResults_10percentileSE; % average - 10th p)
% 
% % Calculate the number of groups and number of bars in each group
% [ngroups,nbars] = size(AvAgeResultsSE);
% % Get the x coordinate of the bars
% x = nan(nbars, ngroups);
% for i = 1:nbars
%     x(i,:) = b(i).XEndPoints;
% end
% 
% er = errorbar(x', AvAgeResultsSE.', p_errlow.', p_errhigh.', 'k', 'LineWidth',1,'linestyle','none');    

labels = {'Past','1.5','2.0','3.0'};
legend(labels);
% set(legend,...
%     'Position',[0.594364046168392 0.717496545895171 0.116719241821051 0.13589211247274]);
% 
% % Create ylabel
% ylabel('Average Annual Heat Related Deaths per 100,000 population');
% 
% % Create xlabel
% xlabel('Age Group');
% %create title
% title ('Climate and Population Change');
% 
% set (gca, 'ylim', [0 M])
% set(gca,'XTick',[1 2 3 4],'XTickLabel',{'0-64', '65-74', '75-84','85+'},'YGrid',...
%     'on');

%% Clustered bars per age group, climate scenario, cc and cc+SE
% Absolute deaths

% figure
% M = max((AvAnnualDeathsUK_SE_age_90th(:,4))*1.1);
% subplot(1,2,1);
% b = bar(AvAnnualDeathsUK_age.'); %create bar
% hold on
% 
% set (gca, 'ylim', [0 M])
% set(gca, 'YTickLabel',get(gca,'YTick')) % avoid exponential on legend
% set(gca,'XTick',[1 2 3 4],'XTickLabel',{'0-64', '65-74', '75-84','85+'},'YGrid',...
%     'on');
% 
% hold on
% p_errhigh = AvAnnualDeathsUK_age_90th - AvAnnualDeathsUK_age; %90th p - average
% p_errlow = AvAnnualDeathsUK_age - AvAnnualDeathsUK_age_10th; % average - 10th p)
% 
% % Calculate the number of groups and number of bars in each group
% [ngroups,nbars] = size(AvAnnualDeathsUK_age);
% % Get the x coordinate of the bars
% x = nan(nbars, ngroups);
% for i = 1:nbars
%     x(i,:) = b(i).XEndPoints;
% end
% 
% er = errorbar(x', AvAnnualDeathsUK_age.', p_errlow.', p_errhigh.', 'k', 'LineWidth',1,'linestyle','none');                                
% 
% % Create ylabel
% ylabel('Average Annual Heat Related Deaths');
% 
% % Create xlabel
% xlabel('Age Group');
% %create title
% title ('Climate Change Only')
% 
% subplot(1,2,2);
% b= bar(AvAnnualDeathsUK_SE_age.'); %create bar
% hold on
% 
% set (gca, 'ylim', [0 M])
% set(gca, 'YTickLabel',get(gca,'YTick')) % avoid exponential on legend
% set(gca,'XTick',[1 2 3 4],'XTickLabel',{'0-64', '65-74', '75-84','85+'},'YGrid',...
%     'on');
% 
% 
% hold on
% p_errhigh = AvAnnualDeathsUK_SE_age_90th - AvAnnualDeathsUK_SE_age; %90th p - average
% p_errlow = AvAnnualDeathsUK_SE_age - AvAnnualDeathsUK_SE_age_10th; % average - 10th p)
% 
% % Calculate the number of groups and number of bars in each group
% [ngroups,nbars] = size(AvAnnualDeathsUK_SE_age);
% % Get the x coordinate of the bars
% x = nan(nbars, ngroups);
% for i = 1:nbars
%     x(i,:) = b(i).XEndPoints;
% end
% 
% er = errorbar(x', AvAnnualDeathsUK_SE_age.', p_errlow.', p_errhigh.', 'k', 'LineWidth',1,'linestyle','none');    
% 
% labels = {'Past','1.5','2.0','3.0'};
% legend(labels);
% set(legend,...
%     'Position',[0.594364046168392 0.717496545895171 0.116719241821051 0.13589211247274]);
% 
% % Create ylabel
% ylabel('Average Annual Heat Related Deaths');
% 
% % Create xlabel
% xlabel('Age Group');
% %create title
% title ('Climate and Population Change');

%% clustered bar - mortality per each 1 degree increment for each scenario
figure
subplot(2,1,1);

avMortPerIncrementUK_0dp = round(avMortPerIncrementUK); % round to nearest integer
M = max(avMortPerIncrementUK_0dp(8,:)*1.1); %set y axis max limit across sub-plots

b = bar(avMortPerIncrementUK_0dp([2:3, 5, 7],1:15).');
set(b(1),'FaceColor',[1 0.843137264251709 0]);
set(b(2),'FaceColor',[0.87058824300766 0.490196079015732 0]);
set(b(3),'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(b(4),'FaceColor',[0.600000023841858 0.200000002980232 0]);

% Create ylabel
ylabel('Average Annual Heat Related Deaths');

% Create xlabel
xlabel('Degrees warming above regional thresholds');
%create title
title ('Climate Change Only')
set (gca, 'ylim', [0 M])


subplot(2,1,2);
b = bar(avMortPerIncrementUK_0dp([2,4, 6, 8],1:15).');
set(b(1),'FaceColor',[1 0.843137264251709 0]);
set(b(2),'FaceColor',[0.87058824300766 0.490196079015732 0]);
set(b(3),'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(b(4),'FaceColor',[0.600000023841858 0.200000002980232 0]);
hold on
labels = {'Past','1.5','2.0','3.0'};
legend(labels);

% Create ylabel
ylabel('Average Annual Heat Related Deaths');
xlabel('Degrees warming above regional thresholds');
%create title
title ('Climate and Population Change')
set (gca, 'ylim', [0 M])

%% Plots for Adaptation (scenarios for natural acclimatisation)
figure;
subplot(1,2,1);
b= bar (AvAnnualDeathsAdapt.', 'stacked', 'FaceColor','flat');
b(1).CData = [0.800000011920929 0.800000011920929 0.800000011920929];
b(2).CData = [0.501960813999176 0.501960813999176 0.501960813999176];
b(3).CData = [0.313725501298904 0.313725501298904 0.313725501298904];
hold on

M = AvAnnualDeathsUK_SE(1,10)*1.1; %set y axis max limit across sub-plots
% Create ylabel
ylabel('Average Annual Heat Related Deaths');

% Create xlabel
xlabel('Climate Scenario');

% Create title
title('Climate Change Only');

 labels = {'Adaptation 2°C','Adaptation 1°C','No Adaptation'};
 legend(labels);
 set(legend,...
     'Position',[0.157187955560562 0.762734641133266 0.276785708750997 0.13589211247274]);

% Set the remaining axes properties
set (gca, 'ylim', [0 M])
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Past','1.5','2.0','3.0'},'YGrid',...
    'on');

subplot(1,2,2);
b=bar (AvAnnualDeathsAdapt_SE.', 'stacked', 'FaceColor','flat');
b(1).CData = [0.800000011920929 0.800000011920929 0.800000011920929];
b(2).CData = [0.501960813999176 0.501960813999176 0.501960813999176];
b(3).CData = [0.313725501298904 0.313725501298904 0.313725501298904];
hold on

% Create ylabel
ylabel('Average Annual Heat Related Deaths');

% Create xlabel
xlabel('Climate Scenario');

% Create title
title('Climate and Population Change');

% Set the remaining axes properties
set (gca, 'ylim', [0 M])
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Past','1.5','2.0','3.0'},'YGrid',...
    'on');

% %% COMPARE bias vs non-bias -TEST ONLY -REMOVE
% M = 15000; %set y axis max limit across sub-plots
% figure
% subplot(1,2,1);
% b = bar(unnamed(1:4, 1:2), 'FaceColor', 'flat');
% b(1).CData(1,:) = [0.494117647409439 0.494117647409439 0.494117647409439];
% b(2).CData(1,:) = [0.235294118523598 0.235294118523598 0.235294118523598];
% b(1).CData(2,:) = [0.494117647409439 0.494117647409439 0.494117647409439];
% b(2).CData(2,:) = [0.235294118523598 0.235294118523598 0.235294118523598];
% b(1).CData(3,:) = [0.494117647409439 0.494117647409439 0.494117647409439];
% b(2).CData(3,:) = [0.235294118523598 0.235294118523598 0.235294118523598];
% b(1).CData(4,:) = [0.494117647409439 0.494117647409439 0.494117647409439];
% b(2).CData(4,:) = [0.235294118523598 0.235294118523598 0.235294118523598];
% 
% hold on
% 
% % Calculate the number of groups and number of bars in each group
% [ngroups,nbars] = size(unnamed(1:4, 1:2));
% % Get the x coordinate of the bars
% x = nan(nbars, ngroups);
% for i = 1:nbars
%     x(i,:) = b(i).XEndPoints;
% end
% 
% er = errorbar(x', unnamed(1:4, 1:2), unnamed(5:8, 1:2), unnamed(9:12, 1:2), 'k', 'LineWidth',1,'linestyle','none');                                
% 
% % Create ylabel
% ylabel('Average Annual Heat Related Deaths');
% 
% % Create xlabel
% xlabel('Climate Scenario');
% 
% % Create title
% title('Climate Change Only');
% 
% % Set the remaining axes properties
% set (gca, 'ylim', [0 M])
% set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Past','1.5','2.0','3.0'},'YGrid',...
%     'on');
% 
% subplot(1,2,2);
% b = bar(unnamed(1:4, 3:4), 'FaceColor', 'flat');
% b(1).CData(1,:) = [0.494117647409439 0.494117647409439 0.494117647409439];
% b(2).CData(1,:) = [0.235294118523598 0.235294118523598 0.235294118523598];
% b(1).CData(2,:) = [0.494117647409439 0.494117647409439 0.494117647409439];
% b(2).CData(2,:) = [0.235294118523598 0.235294118523598 0.235294118523598];
% b(1).CData(3,:) = [0.494117647409439 0.494117647409439 0.494117647409439];
% b(2).CData(3,:) = [0.235294118523598 0.235294118523598 0.235294118523598];
% b(1).CData(4,:) = [0.494117647409439 0.494117647409439 0.494117647409439];
% b(2).CData(4,:) = [0.235294118523598 0.235294118523598 0.235294118523598];
% 
% hold on
% 
% % Calculate the number of groups and number of bars in each group
% [ngroups,nbars] = size(unnamed(1:4, 1:2));
% % Get the x coordinate of the bars
% x = nan(nbars, ngroups);
% for i = 1:nbars
%     x(i,:) = b(i).XEndPoints;
% end
% 
% er = errorbar(x', unnamed(1:4, 3:4), unnamed(5:8, 3:4), unnamed(9:12, 3:4), 'k', 'LineWidth',1,'linestyle','none');            
% % Create ylabel
% ylabel('Average Annual Heat Related Deaths');
% 
% % Create xlabel
% xlabel('Climate Scenario');
% 
% % Create title
% title('Climate and Population Change');
% 
% % Set the remaining axes properties
% set (gca, 'ylim', [0 M])
% set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Past','1.5','2.0','3.0'},'YGrid',...
%     'on');
elseif impactMetric == 2
    run ARCADIA_HEAT_Impacts_labour_productivity.m
end
toc  %stopwatch to measure performance
disp END