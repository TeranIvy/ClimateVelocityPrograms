% % % ==================================================================== % % %
% % % 
% % % savedata_serial.m
% % % 
% % % Author: Iva Kavcic
% % %
% % % Date last modified: 17/11/2016
% % % Date of last comments update: 08/04/2018
% % % Runs with Matlab versions R2010a and newer
% % % 
% % % This program reads the monthly maximum temperature, minimum
% % % temperature and total precipitation data in ASCII format and
% % % saves them as binary data stored in *.mat files. It can save
% % % yearly or seasonal data.
% % % The data can be downloaded from
% % % ftp://ftp-cdc.dwd.de/pub/CDC/grids_germany/monthly/
% % % and for terms of use please see
% % % ftp://ftp-cdc.dwd.de/pub/CDC/Terms_of_use.pdf.
% % % The original data are stored in integer format, 1/10 degree
% % % Celsius for temperature and mm for precipitation (e.g.
% % % -57 °C is -5.7 °C).
% % % Instead of dividing the temperature data by 10 before storing them,
% % % their integer format is preserved in order to be able to use single
% % % precision for the generated binary storage files and save space.
% % % Single precision gives 6 reliable decimal digits so the quality
% % % of stored data is not degraded.
% % % 
% % % Additional functionality of this program is to calculate data
% % % averages on 10 km by 10 km grid and store them to output 
% % % binary *.mat files. The averaged data are also stored in single
% % % precision for consistency (meteorological data are usually
% % % stored with 1 or 2 decimal digits). 
% % %  
% % % ==================================================================== % % %
%
% ---------- Define required parameters to process the data files ------------ %
%
% Clear workspace
close all; clear all;
%
% Number of parameters in input ASCII file which describe the data
nvar = 6;
% Original data are stored in (NROWS, NCOLS) matrix in monthly ASCII
% input data files. Since values of NROWS and NCOLS are same for all
%  files, they are hardcoded here
NROWS = 866;
NCOLS = 654;
%
% Averaging step (number of data in each direction)
naverage = 10;
navstr = num2str(naverage); 
% Calculate size of averaged data
NROWS_av = floor(NROWS/naverage);
NCOLS_av = floor(NCOLS/naverage);
nrowdat = NROWS_av*naverage;
ncolsdat = NCOLS_av*naverage;
NROWS_rest = NROWS - nrowdat;
NCOLS_rest = NCOLS - ncolsdat;
% Decide which data to cut from the averaging process - here even number 
% of columns from the left and right and rows from top and bottom.
% If there is an odd number of rows and columns to cut, one more row will
% be cut from top and one more column from the left.
% As seen from inspecting ASCII datafiles, these rows and columns are 
% usually filled with -999 (missing data) so the loss is not significant.
NROWS_up = ceil(NROWS_rest/2);     NROWS_down = NROWS_rest - NROWS_up;
NCOLS_left = ceil(NCOLS_rest/2);   NCOLS_right = NCOLS_rest - NCOLS_left;
%
% Home directory path
pathhome = '..\';
% Temperature and precipitation data are stored in separate
% directories (vardir).
vardir = {'MaxTemp','MinTemp','Precipitation'};
% Capture variable descriptors of the input ASCII data files: 
% TADXM for maximum temperature, TADNMM for minimum temperature
% and RSMS for precipitation. For instance, the data files containing
% information for February 1910 are stored in TADXMM_02_1910_01.asc,
% TADNMM_02_1910_01.asc and  RSMS_02_1910_01.asc (the last "01" string
% denotes resolution of data).
varfile = {'TADXMM','TADNMM','RSMS'};
%
% Define range of years available to process. Here the data range available
% is January 1901 to February 2016 (Winter season)
yearmin = 1901;   yearmax = 2015;
%
% List months in a year and define seasons
months_all = {'01','02','03','04','05','06','07','08','09','10','11','12'};
% Months in all seasons (seasons are defined following climate research
% conventions, e.g. Dec-Jan-Feb for Winter)
months_seas{1} = {months_all(12),months_all(1),months_all(2)};   % Winter
months_seas{2} = months_all(3:5);                                % Spring
months_seas{3} = months_all(6:8);                                % Summer
months_seas{4} = months_all(9:11);                               % Autumn
names_seas = {'Win','Spr','Sum','Aut'};
%
% ---------- Select variable, years and season to process  ------------------- %
%
% Choose the variable to process
display('Choose the variable to process: 1 for Maximum Temperature,');
display('                                2 for Minimum Temperature or');
ivarb = input('                                3 for Precipitation: ');
display(' ');
% Set input and output paths
pathout = [pathhome 'Grids_Germany_' char(vardir(ivarb)) '_GZ\'];
pathin = [pathout 'Files_ASC_' char(vardir(ivarb)) '\'];
%
% Choose the start and end year to process
yearstart = input('Choose the start year to process (1901 - 2015): ');
yearend   = input('Choose the end year to process (1901 - 2015): ');
%
% Choose the season to process
display('Choose the season to process: 1 for Winter,');
display('                              2 for Spring,');
display('                              3 for Summer,');
display('                              4 for Autumn or');
iseas = input('                              any other number for the entire year: ');
display(' ');
%
% Make range of years and months to process
years_index = [yearstart:yearend];   nyears = length(years_index);
yyyy_start = num2str(years_index(1));
yyyy_end = num2str(years_index(end));
% All months for processing yearly data or months selected from
% months_seas if processing seasonall data
if iseas >= 1 && iseas <= 4;
   months_str = months_seas{iseas};
   fseas = char(names_seas{iseas});
else;
   months_str = months_all;
   fseas = 'All';
end;
nmonths = length(months_str);
%
% Allocate some memory for the data (single precision for output
% binary data files)
ndata = nyears*nmonths;   % Length of time series for each grid point
% Individual 1 km data are output to (NROWS,NCOLS,ndata) matrix
data1km_all = single(zeros(NROWS,NCOLS,ndata));
% Averaged 10 km data are output to (NROWS_av,NCOLS_av,ndata) matrix
data10km_all = single(zeros(NROWS_av,NCOLS_av,ndata));
% Dates are output to (year,month) array of length ndata
dates_all = single(zeros(ndata,2));
%
% ---------- Process data on a year-by-year basis  --------------------------- %
%
% Loop over the selected years
for iy = 1:nyears;
    % Determine current and next year
    yearcurr = years_index(iy);   yearnext = yearcurr + 1;
    yyyy = num2str(yearcurr);     yyyp = num2str(yearnext);
    ystr = repmat({yyyy},1,nmonths);
    if iseas == 1;
       % For Winter we need to process January and February of the next year
       ystr(end-1:end) = {yyyp};
    end;
    % Loop over the selected months in a year
    for im = 1:nmonths;
        % Make input filename (e.g. RSMS_03_1950_01.asc)
        km = (iy - 1)*nmonths + im;
        mm = char(months_str{im});   yy = char(ystr{im});
        dates_all(km,1) = str2num(mm);
        dates_all(km,2) = str2num(yy);
        fnamein = [pathin char(varfile(ivarb)),'_' mm '_' yy '_01.asc'];
        % Open input ASCII data file with 1 km values
        fidin = fopen(fnamein,'r');
        % Read first 6 (nvar) lines with control variables and convert them
        % to numerical values
        variable_names = cell(1,nvar);
        variable_values = cell(1,nvar);
        for i = 1:nvar;
            hline = fgetl(fidin);
            varstr = textscan(hline,'%s %s');
            variable_names{i} = char(varstr{1});
            variable_values{i} = char(varstr{2});    
            v = genvarname(variable_names{i});
            eval([v ' = str2num(variable_values{i});']);
        end;
        % Data are stored in (NROWS,NCOLS) matrix after the line where
        % control variables are listed
        % Initialise temporary data matrix to 0
        data1km = single(zeros(NROWS,NCOLS));
        % Read data row by row and convert them to numerical values
        for irow = 1:NROWS;
            tline = fgetl(fidin);
            data1km(irow,:) = single(str2num(tline));
        end;
        % Close input ASCII data file
        fclose(fidin);
        % Replace missing data with NaN (NODATA_VALUE is -999)
        data1km(find(data1km == NODATA_VALUE)) = NaN;
        % Matlab reads the northernmost row as the first (lowest index)
        % so the temporary data field needs to be flipped back to the
        % correct orientation as in the input file (southernmost row
        % having the lowest index)
        data1km = flipud(data1km);
        % Create time index (1 to ndata) for one month of data
        kym = (iy - 1)*nmonths + im;
        % Store 1 km values
        data1km_all(:,:,kym) = data1km;
        % Extract only the data to average
        datamat = data1km(NROWS_up+1:end-NROWS_down,...
                          NCOLS_left+1:end-NCOLS_right);
        % Calculate spatial (10 km)^2 averages using Matlab function mean2
        data10km = single(NaN*ones(NROWS_av,NCOLS_av));
        for irow = 1:NROWS_av;
            i1 = (irow - 1)*naverage + 1;   i2 = irow*naverage;
            datastrip = datamat(i1:i2,:);
            for jcol = 1:NCOLS_av;
                j1 = (jcol - 1)*naverage + 1;   j2 = jcol*naverage;
                data10km(irow,jcol) = single(mean2(datastrip(:,j1:j2)));
            end;
        end;
        % Store values of averaged data
        data10km_all(:,:,kym) = data10km;
        % Clear some temporary inner loop variables
        clear datamat i irow jcol i1 i2 j1 j2 datastrip fidin fnamein; 
    end;
    % Display which year was processed
    display(['Processed year ' yyyy]);
    % Clear some temporary outer loop variables
    clear yearcurr yearnext yyyy yyyp ystr;
end;
%
% Calculate cellsize (cell length) for averaged data and store it. 
% Here the cell length of averaged data is 10 km and is the same
% in each direction.
CELLSIZE10 = CELLSIZE*naverage;
%
% ---------- Save data in binary *.mat format  ------------------------------- %
%
% Create output file names to resemble the format of the input ASCII file
% file names. For example, the binary mat file containing values of
% maximum temperatures for Spring from 1901 to 1950 will be called
% 'TADXMM_1901_1950_Spr.mat' and the *.mat file containing values for an
% entire year from 1951 to 1975 will be called 'TADXMM_1951_1975_All.mat'.
fstr = [char(varfile(ivarb)) '_' yyyy_start '_' yyyy_end];
%
% Save 1 km original data in *.mat format
fnameout_1km = [pathout fstr '_' fseas '.mat'];
save(fnameout_1km, 'data1km_all', 'dates_all', 'NROWS', 'NCOLS', ... 
    'XLLCORNER', 'YLLCORNER', 'NODATA_VALUE', 'CELLSIZE', '-single');
%
% Save 10 km averaged data in *.mat format (additional '_av10' in the file
% name, e.g. 'TADXMM_1901_1950_av10_Spr.mat'
fnameout_10km = [pathout fstr '_av' navstr '_' fseas '.mat'];
save(fnameout_10km, 'data10km_all', 'dates_all', 'NROWS_av', 'NCOLS_av', ...
    'XLLCORNER', 'YLLCORNER', 'NODATA_VALUE', 'CELLSIZE10', '-single');
	