% % % ==================================================================== % % %
% % % 
% % % SpatialGradients_Season.m
% % % 
% % % Author: Iva Kavcic
% % % 
% % % Date last modified: 19/05/2017
% % % Date of last comments update: 08/04/2018
% % % Runs with Matlab versions R2010a and newer
% % % 
% % % This program reads in the 1 km data stored in binary form in *.mat
% % % files and then calculates the median spatial gradients by calling
% % % function SpatGrad_Median.m.
% % %
% % % Input *.mat data files containing 1 km seasonal data were created
% % % using savedata_serial.m script for each analysed variable (maximum
% % % temperature, minimum temperature and precipitation). Two files
% % % were created for each variable and season, covering periods
% % % [1901 1950] and [1951 2015].
% % % Spatial gradients can be calculated for three time periods:
% % % [1901 1950], [1951 2015] and [1901 2015]. 
% % %
% % % The results are initialised to NaN and stored in binary *.mat files
% % % in double precision.
% % %  
% % % ==================================================================== % % %
%
% ---------- Define required parameters to process the data files ------------ %
%
% Clear workspace
close all; clear all;
%
% Home directory path
pathhome = '..\';
%
% Temperature and precipitation data are stored in separate
% directories (vardir).
vardir = {'MaxTemp','MinTemp','Precipitation'};
% Capture variable descriptors of the input *.mat data files: 
% TADXM for maximum temperature, TADNMM for minimum temperature and RSMS
% for precipitation. For example, the data file containing values of maximum
% temperatures for Spring from 1951 to 2015 is called
% 'TADXMM_1951_2015_Spr.mat'.
varfile = {'TADXMM','TADNMM','RSMS'};
names_seas = {'Win','Spr','Sum','Aut'};
% The original data are stored in integer format, 1/10 degree Celsius
% for temperature and mm for precipitation, so temperature data need to
% be divided by 10. To avoid "if" test, scaling factors are defined here
% in an array so the one for the selected variable can be used later.
vardivscale = [10, 10, 1];
%
% Set identifiers for output paths. The calculated spatial gradients values
% will be stored in relevant subdirectories within the "pathin" directories.
dirout = 'Spatial_Gradients/';
% Set the identifier for output file names
spgrtext = 'SpatialGradient';
%
% Define range of years available to process according to the time
% periods of the input *.mat data files
years_process = [1901 1950; 1951 2015];
%
% ---------- Select variable, period and season to process ------------------- %
% (If non-interactive run is required, ivarb, iper and iseas can be hardcoded)
%
% Choose the variable to process
display('Choose the variable to process: 1 for Maximum Temperature,');
display('                                2 for Minimum Temperature or');
ivarb = input('                                3 for Precipitation: ');
display(' ');
% Set input and output paths
pathin = [pathhome 'Grids_Germany_' char(vardir(ivarb)) '_GZ/'];
pathout = [pathin dirout];
% Set portion of the input data file name
varproc = char(varfile(ivarb));
% Set scaling factor
varscale = vardivscale(ivarb);
%
% Choose the period to process
display('Choose the period to process: 1 for control (1901 - 1950),');
display('                              2 for recent (1951 - 2015) or');
iper = input('                              any other number for the entire period: ');
display(' ');
%
% Choose the season to process
display('Choose the season to process: 1 for Winter,');
display('                              2 for Spring,');
display('                              3 for Summer or ');
iseas = input('                              4 for Autumn: ');
display(' ');
% Determine the input season string
if iseas >= 1 && iseas <= 4;
   fseas = char(names_seas{iseas});
else;
   display(' Wrong iseas parameter!');
end;
%
% ---------- Calculating spatial gradients for the whole domain -------------- %
%
% Display information
display(' ')
display(['**********   Spatial gradients for 1 km seasonal data   **********']);
display(' ')
%
% Load the files with all data
if (iper == 1 || iper == 2);
    % For [1901 1950] and [1951 2015] periods load one data file
    yyyy_start = num2str(years_process(iper,1));
    yyyy_end = num2str(years_process(iper,2));
    fstr = [varproc '_' yyyy_start '_' yyyy_end];
    fnamein_1km = [pathin fstr '_' fseas];
    eval(['load ' fnamein_1km ';']);
    % Convert input data (single precision) to double precision
    % before further analysis
    data1km_total = double(data1km_all);
    dates_total = double(dates_all);
    % Clear input data from memory
    clear data1km_all dates_all;
else
    % For [1901 2015] period load both data files
    % Create work arrays
    data1km_total = [];
    dates_total = [];
    for i = 1:2;
        yyyy_start = num2str(years_process(i,1));
        yyyy_end = num2str(years_process(i,2));
        fstr = [varproc '_' yyyy_start '_' yyyy_end];
        fnamein_1km = [pathin fstr '_' fseas];
        eval(['load ' fnamein_1km ';']);
        % Convert input data (single precision) to double precision
        % before further analysis
        data1km_total = cat(3, data1km_total, double(data1km_all)); 
        dates_total = [dates_total; double(dates_all)];
        % Clear input data from memory
        clear data1km_all dates_all;
    end;
    % Redefine start and end years for output files
    yyyy_start = num2str(years_process(1,1));
    yyyy_end = num2str(years_process(end,end));
end;
%
% Scale data with the appropriate factor
data1km_total = data1km_total./varscale;
%
% Get matrix dimensions (number of rows and columns)
[nrows1, ncols1] = size(data1km_total(:,:,1));
%
% Calculate year coordinate as year + (month-1)/12
MyYears = dates_total(:,2) + (dates_total(:,1)-1)./12;
% 
% Find time series of spatial gradients of data 
[sgradmed1km_NS, sgradmed1km_WE, ...
        sgradmed1km_abs, sgradmed1km_angle, inddata1km] = ...
        SpatGrad_Median(data1km_total, CELLSIZE); 
display('Calculated medians of time series of spatial gradients');
%
% ---------- Save spatial gradients in binary *.mat format ------------------- %
%
% Create output file names to resemble the format of the input *.mat files
% For example, the binary *.mat file containing spatial gradients for Summer
% minimum temperature from 1951 to 2015 will be called 
% 'SpatialGradient_TADNMM_1951_2015_Sum.mat'.
fstr = [spgrtext '_' varproc '_' yyyy_start '_' yyyy_end];
fnameout_1km = [dirout fstr '_' fseas '.mat'];
save(fnameout_1km, 'sgradmed1km_NS', 'sgradmed1km_WE', 'sgradmed1km_abs', ...
                   'sgradmed1km_angle', 'inddata1km', 'CELLSIZE');
display('Saved spatial gradients outputs for 1 km seasonal data');
