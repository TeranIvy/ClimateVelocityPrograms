% % % ==================================================================== % % %
% % % 
% % % MKyear_ktaub.m
% % % 
% % % Author: Iva Kavcic
% % % 
% % % Date last modified: 17/11/2016
% % % Date of last comments update: 08/04/2018
% % % Runs with Matlab versions R2010a and newer
% % % 
% % % This program reads in the 1 km data stored in binary form in *.mat
% % % files and then performs Mann-Kendall (hereafter denoted as MK) trend
% % % analysis using Jeff Burkey's Matlab script ktaub.m.
% % % The ktaub.m script, together with the test data and references
% % % describing the methodology, can be retrieved from
% % % http://www.mathworks.com/matlabcentral/fileexchange/authors/23983
% % %
% % % Input *.mat data files containing 1 km yearly data were created
% % % using savedata_serial.m script for each analysed variable (maximum
% % % temperature, minimum temperature and precipitation). Four files
% % % were created for each variable to make them manageable for
% % % analysis, covering periods [1901 1925], [1926 1950], [1951 1975]
% % % and [1976 2015].
% % % MK analysis can be performed for three time periods: [1901 1950],
% % % [1951 2015] and [1901 2015].
% % %
% % % ktaub function performs MK analysis for one time series. Here we
% % % have NROWS*NCOLS time series. Due to computational costs, the
% % % (NROWS,NCOLS,ndata) data matrix is processed in 100-rows strips
% % % for each complete time series (containing all ndata values).
% % %
% % % The data for the selected period are are analysed three times using
% % % different values of significance levels alpha = 0.05, 0.01 and 0.001.
% % % First analysis is performed for all complete time series at the
% % % significance level alpha = 0.05. After that the ktaub return
% % % variable "sig" (p value) is used to determine if a specific time
% % % series will be analysed again for the next significance level
% % % (sig <= alpha) or omitted (sig > alpha).
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
% We do not want plot in MK testing so set ktaub "wantplot" parameter to 0
wantplot = 0;
% 
% Define confidence levels for MK testing (in %)
confidences_all = [95 99 99.9]';
confstr = {'950','990','999'};
nconf = length(confidences_all);
% Calculate significance levels from confidence levels
alphas_all = 1 - confidences_all./100;
%
% Home directory path
pathhome = '..\';
% Temperature and precipitation data are stored in separate
% directories (vardir).
vardir = {'MaxTemp','MinTemp','Precipitation'};
% Capture variable descriptors of the input *.mat data files: 
% TADXM for maximum temperature, TADNMM for minimum temperature and RSMS
% for precipitation. For example, the data file containing values of maximum
% temperatures for the entire year from 1951 to 1975 is called
% 'TADXMM_1951_1975_All.mat'.
varfile = {'TADXMM','TADNMM','RSMS'};
fseas = 'All';
% The original data are stored in integer format, 1/10 degree Celsius
% for temperature and mm for precipitation, so temperature data need to
% be divided by 10. To avoid "if" test, scaling factors are defined here
% in an array so the one for the selected variable can be used later.
vardivscale = [10, 10, 1];
%
% Set identifiers for output paths. The calculated MK test values will be
% stored in relevant subdirectories within the "pathin" directories.
dirout = 'MKparams_ktaub/';
% Set the identifier for output file names
mktext = ['MKktaub'];
%
% Define range of years available to process according to the time
% periods of the input *.mat data files
years_process = [1901 1925; 1926 1950; 1951 1975; 1976 2015];
iper_skip = [0 2];
%
% ---------- Select variable and period to process --------------------------- %
% (If non-interactive run is required, ivarb and iper can be hardcoded)
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
% ---------- MK analysis on a divided domain (100-rows strips) --------------- %
%
% Display information
display(' ')
display('**********   MK test for 1 km year data (ktaub)   **********');
display(' ')
%
% Load first of the data files to get matrix dimensions
% Determine years identifiers for the input data file
yyyy_start = num2str(years_process(1,1));
yyyy_end = num2str(years_process(1,2));
% Create input file name
fstr = [char(varfile(ivarb)) '_' yyyy_start '_' yyyy_end];
fnamein_1km = [pathin fstr '_' fseas];
% Load input file and get matrix dimensions (number of rows and columns)
eval(['load ' fnamein_1km ';']);
[nrows1, ncols1] = size(data1km_all(:,:,1));
% Clear data
clear data1km_all;
%
% Create start and end row indices of strips to process using matrix
% dimensions (nrows1,ncols1)
nrowslice = 100;
nslice = floor(nrows1/nrowslice);
nrowfloor = nslice*nrowslice;
vecslice = [[0:nrowslice:nrowfloor] nrows1];
% Number of slices
nvsl = length(vecslice);
% Create start and end indices for strips
nstrips = nvsl-1;
istrip = zeros(nstrips,2);
for k = 1:nstrips;
    istrip(k,1) = vecslice(k) + 1;
    istrip(k,2) = vecslice(k+1);
end;
%
% Create initialisation matrix filled with NaNs
nan1km = NaN.*ones(nrows1,ncols1);
%
% Now allocate matrices for results of MK analysis using matrix dimensions
% (nrows1,ncols1) and set them to NaN.
% ktaub return parameters - only 3 outputs change with different
% significance levels: h, CIlower and CIupper. They are preallocated
% separately and stored together with the rest of results.
% C3 is an array used to generate Sen slopes and senplot can be calculated
% later in the plotting routine of ktaub.m. Hence, no memory is allocated
% for them here.
%
% Initialisation for alpha = 0.05 (confidence level = 95%)
taub1 = nan1km;          tau1 = nan1km;            h1_950 = nan1km;
sig1 = nan1km;           Z1 = nan1km;              S1 = nan1km;
sigma1 = nan1km;         sen1 = nan1km;            n1 = nan1km; 
CIlower1_950 = nan1km;   CIupper1_950 = nan1km;    D1 = nan1km;
Dall1 = nan1km;          nsigma1 = nan1km;
% Initialisation for alpha = 0.01 (confidence level = 99%)
h1_990 = nan1km;         CIlower1_990 = nan1km;    CIupper1_990 = nan1km; 
% Initialisation for alpha = 0.001 (confidence level = 99.9%)
h1_999 = nan1km;         CIlower1_999 = nan1km;    CIupper1_999 = nan1km; 
%
% Load the files with data and do MK test for each strip in turn
for kstrip = 1:nstrips;
    % Start and end row
    istart = istrip(kstrip,1);   iend = istrip(kstrip,2);
    % Create work arrays
    data1km_total = [];
    dates_total = [];
    % Load the files with data
    if (iper == 1 || iper == 2);
        % For [1901 1950] and [1951 2015] periods load two data files
        for i = 1:2;
            k = iper_skip(iper) + i;
            yyyy_start = num2str(years_process(iper,1));
            yyyy_end = num2str(years_process(iper,2));
            fstr = [varproc '_' yyyy_start '_' yyyy_end];
            fnamein_1km = [pathin fstr '_' fseas];
            eval(['load ' fnamein_1km ';']);
            % Convert input data (single precision) to double precision
            % before further analysis
            data1km_total = cat(3, data1km_total, ...
                double(data1km_all(istart:iend,:,:)));
            dates_total = [dates_total; double(dates_all)];
            % Clear input data from memory
            clear data1km_all dates_all; 
        end;         
    else;
        % For [1901 2015] period load all four data files
        for i = 1:4;
	        yyyy_start = num2str(years_process(i,1));
            yyyy_end = num2str(years_process(i,2));
            fstr = [varproc '_' yyyy_start '_' yyyy_end];
            fnamein_1km = [pathout fstr '_' fseas];
            eval(['load ' fnamein_1km ';']);
            % Convert input data (single precision) to double precision
            % before further analysis
            data1km_total = cat(3, data1km_total, ...
                double(data1km_all(istart:iend,:,:)));
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
    % Length of time series
    ndata = size(data1km_total,3);
    % Calculate year coordinate as year + (month-1)/12
    MyYears = dates_total(:,2) + (dates_total(:,1)-1)./12;
    display(['Loaded 1 km data for ' num2str(kstrip) ...
	         '. strip, now doing MK test: ']);
    %
    % Call MK test with wantplot = 0 and alpha = 0.05 (confidence level 95%)
    % for all points in this strip of data1km_total
    alpha_sig = alphas_all(1);
    display(['MK ktaub test for confidence = ' ...
              num2str(confidences_all(1)) '% : ']);
    % Initialise work arrays	
    notnans1 = nan1km;
    a = zeros(ndata,1);
    % Work row by row in the strip
    for i = istart:iend;
        iloc = i - istart + 1;
        for j = 1:ncols1;
            % Isolate one time series and check whether it is complete
            % (not containing NaNs)
            a(:) = data1km_total(iloc,j,:);
            notnans1(i,j) = ~any(isnan(a));
            % If the time series is complete call ktaub function
            if (notnans1(i,j) == 1);
               datain = [MyYears a];
               [taub1(i,j) tau1(i,j) h1_950(i,j) sig1(i,j) Z1(i,j) S1(i,j)...
                sigma1(i,j) sen1(i,j) n1(i,j) senplot ...
                CIlower1_950(i,j) CIupper1_950(i,j) D1(i,j) ...
                Dall1(i,j) C3 nsigma1(i,j)] = ...
                ktaub(datain, alpha_sig, wantplot);
            end;
        end;
        display([' Processed row ', num2str(i),', iloc = ',num2str(iloc)]);
    end;
    display(' ');
    %
    % Call MK test with wantplot = 0 and alpha = 0.01 (confidence level 99%)
    alpha_sig = alphas_all(2);
    % Initialise work array
    a = zeros(ndata,1);
    % Find where sig1 <= 0.01 and analyse again only those time series
    % to calculate relevant h, CIlower and CIupper
    [row990, col990] = find(sig1(istart:iend,:) <= alpha_sig);
    if ~isempty(row990);
       nvalues = length(row990);
       display(['MK ktaub test for confidence = ' ...
                 num2str(confidences_all(2)) '% : ']);
       for k = 1:nvalues;
           iloc = row990(k);  j = col990(k);
           i = iloc + istart - 1;
           a(:) = data1km_total(iloc,j,:);
           datain = [MyYears a];
           [taub tau h1_990(i,j) sig Z S sigma sen n senplot ...
            CIlower1_990(i,j) CIupper1_990(i,j) D Dall C3 nsigma] = ...
            ktaub(datain, alpha_sig, wantplot);
        end;
        display(' ');
    end; 
    %
    % Call MK test with wantplot = 0 and alpha = 0.001 (confidence level 99.9%)
    % Do MK ktaub for confidence 99.9%
    alpha_sig = alphas_all(3);
    % Initialise work array
    a = zeros(ndata,1);
    % Find where sig1 <= 0.001 and analyse again only those time series
    % to calculate relevant h, CIlower and CIupper
    [row999, col999] = find(sig1(istart:iend,:) <= alpha_sig); 
    if ~isempty(row999);
       nvalues = length(row999);
       display(['MK ktaub test for confidence = ' ...
                 num2str(confidences_all(3)) '% : ']);
       for k = 1:nvalues;
           iloc = row999(k);  j = col999(k);
           i = iloc + istart - 1;
           a(:) = data1km_total(iloc,j,:);
           datain = [MyYears a];
           [taub tau h1_999(i,j) sig Z S sigma sen n senplot ...
            CIlower1_999(i,j) CIupper1_999(i,j) D Dall C3 nsigma] = ...
            ktaub(datain, alpha_sig, wantplot);
        end;
        display(' ');
    end;
    % Clear some work arrays
    clear row990 col990 row999 col999;
end;
%
% ---------- Save MK parameters in binary *.mat format ----------------------- %
%
% Create output file names to resemble the format of the input *.mat files
% For example, the binary *.mat file containing MK analysis outputs for
% maximum temperature from 1951 to 2015 will be called 
% 'MKktaub_TADXMM_1951_2015_All.mat'.
fstr = [varproc '_' yyyy_start '_' yyyy_end];
fnameout_1km = [pathout mktext '_' fstr '_' fseas '.mat'];
% Save 1 km MK outputs for all analysed significance levels
save(fnameout_1km, 'taub1', 'tau1', 'h1_950', 'sig1', 'Z1', 'S1', 'sigma1', ...
    'sen1', 'n1', 'CIlower1_950', 'CIupper1_950', 'D1', 'Dall1', 'nsigma1', ...
    'h1_990', 'CIlower1_990', 'CIupper1_990', ...
    'h1_999', 'CIlower1_999', 'CIupper1_999', ...
    'confidences_all', 'alphas_all');
display('Saved MK outputs for 1 km year data');
