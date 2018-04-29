% % % ==================================================================== % % %
% % % 
% % % SpatialGradUncertaintyMatlab_Season.m
% % % 
% % % Author: Iva Kavcic
% % % 
% % % Date last modified: 19/05/2017
% % % Date of last comments update: 09/04/2018
% % % Runs with Matlab versions R2010a and newer
% % % 
% % % This program reads in the 1 km data stored in binary form in *.mat
% % % files and then calculates the uncertainties in median spatial
% % % gradients. Before calculating the spatial gradients the data are
% % % detrended using Sen slopes from Mann-Kendall trend analysis.
% % %
% % % Required functions: Detrend_Sen.m        
% % %                     SpatGrad_Median_Unc.m  
% % %                     SpatGradMatlab_Limits_Unc.m
% % %
% % % Input *.mat data files containing 1 km seasonal data were created
% % % using savedata_serial.m script for each analysed variable (maximum
% % % temperature, minimum temperature and precipitation). Two files
% % % were created for each variable and season, covering periods
% % % [1901 1950] and [1951 2015].
% % % Uncertainties in spatial gradients can be calculated for three time
% % % periods: [1901 1950], [1951 2015] and [1901 2015]. 
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
% Set identifiers for input paths to files containing Mann-Kendall results
dirinMK = 'MKparams_ktaubMult/';
mktext = 'MKktaubMult'; 
% Set identifiers for output paths. The calculated uncertainties of spatial
% gradients will be stored in relevant subdirectories within the "pathin"
% directories.
dirout = 'Spatial_Gradients/';
% Set the identifier for output file names
spgrunctext = 'SpatGradUncMatlab';
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
pathinMK = [pathin dirinMK];
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
% --- Calculating spatial gradients with uncertainties for the whole domain -- &
%
% Display information
display(' ')
display(['*****   Spatial gradients for 1 km seasonal detrended data   *****']);
display(['*****************   with uncertainty estimates   *****************']);
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
% ---------- Load file with MK trend analysis results ------------------------ %
%
% Define input path, file name and load the input file
fstr = [mktext '_' varproc '_' yyyy_start '_' yyyy_end]; 
fnameinMK_1km = [pathinMK fstr '_' fseas];
eval(['load ' fnameinMK_1km ';']);
%
% Set significance level for MK results (alphas_all are read from MK file)
alpha_sig = alphas_all(1);
%
% Clear unused fields (we only do alpha = 0.05 here)
clear taub1 tau1 Z1 S1 sigma1 n1 CIlower1_950 CIupper1_950 ...
      D1 Dall1 nsigma1 h1_990 CIlower1_990 CIupper1_990 ...
      h1_999 CIlower1_999 CIupper1_999
%
% ---------- Detrend data and calculate spatial gradients uncertainties ------ %
%
% Detrend data using Sen slope and time coordinate
data1km_detr = Detrend_Sen(data1km_total, sen1, MyYears);
display('Detrended data');
% Free some memory
clear data1km_total
%
% Find time series of spatial gradients on detrended data        
[sgradmedDMB1km_abs, sgradhasdata1km_abs, inddataDMB1km, ...
    nrowssg1km, ncolssg1km] = SpatGrad_Median_Unc(data1km_detr, CELLSIZE); 
display('Calculated median and time series of spatial gradients');
%
% Find limits and spread
[sgradmedDMB1km_limits, sgradmedDMB1km_stdev] = ...
    SpatGradMatlab_Limits_Unc(sgradhasdata1km_abs, nrowssg1km, ...
    ncolssg1km, inddataDMB1km);
display('Found limits and spread with bootstrapping');
%
% ---------- Save spatial gradients uncertainties in binary *.mat format ----- %
%
% Create output file names to resemble the format of the input *.mat files
% For example, the binary *.mat file containing spatial gradients uncertainties
% for Summer minimum temperature from 1951 to 2015 will be called 
% 'SpatGradUncMatlab_TADNMM_1951_2015_Sum.mat'.
fstr = [spgrunctext '_' varproc '_' yyyy_start '_' yyyy_end];
fnameout_1km = [dirout fstr '_' fseas '.mat'];
save(fnameout_1km, 'sgradmedDMB1km_abs', 'sgradmedDMB1km_limits', ...
                   'sgradmedDMB1km_stdev', 'inddataDMB1km');
display('Saved spatial gradients uncertainties outputs for 1 km seasonal detrended data');
