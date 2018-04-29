% % % ==================================================================== % % %
% % % 
% % % function SpatGrad_Median_Unc.m
% % % 
% % % Author: Iva Kavcic
% % % 
% % % Date last modified: 19/05/2017
% % % Date of last comments update: 09/04/2018
% % % Runs with Matlab versions R2010a and newer
% % % 
% % % This function uses 1 km seasonal data read by the main program
% % % SpatialGradUncertaintyMatlab_Season.m to calculate the median and
% % % complete time series of absolute spatial gradients required for
% % % the uncertainty analysis later. The data are first detrended using
% % % function Detrend_Sen.m.
% % %
% % % Input variables:
% % %       datain - Input data
% % %       CELLSIZE - Size of spatial cell
% % % Output variables:   
% % %       sgradmed_abs - Absolute values of spatial gradients, medians
% % %                      over years, dimension [nrowssg, ncolssg]
% % %       sgradhasdata_abs - Reshaped time series of absolute spatial
% % %                          gradients for the uncertainty analysis,
% % %                          dimension [ndata,mdata]
% % %       inddata - Indices of rows and columns for complete time
% % %                 series, dimension [mdata,2]
% % %       nrowssg - Number of rows of 2D median spatial gradient field
% % %       ncolssg - Number of columns of 2D median spatial gradient
% % %                 field
% % %
% % % Spatial gradients are calculated by the 9-points method described in
% % % Burrows et al. 2011, The Pace of Shifting Climate in Marine and
% % % Terrestrial Ecosystems (Science 334:652-655),
% % % http://science.sciencemag.org/content/334/6056/652
% % %  
% % % ==================================================================== % % %
%
function [sgradmed_abs, sgradhasdata_abs, inddata, nrowssg, ncolssg] = ...
          SpatGrad_Median_Unc(datain, CELLSIZE)
    %
    % Calculate size of input data
    [nrows, ncols, ndata] = size(datain);
    % Calculate size of spatial gradients data
    nrowssg = nrows - 2;    ncolssg = ncols - 2;
    %
    % Exclude all non-complete time series
    % First calculate number of data for each grid point
    temp = sum(~isnan(datain),3);
    % Set the control field to NaN if there are missing data
    temp(temp < ndata) = NaN;
    hasalldata = temp(2:end-1,2:end-1);
    % Find row and column indices with the complete time series
    [inddata(:,1), inddata(:,2)] = find(~isnan(hasalldata)); 
    mdata = size(inddata,1);
    %
	% Initialise output fields
    sgradhasdata_abs = NaN.*ones(ndata,mdata);
    sgradmed_abs = NaN.*ones(nrowssg,ncolssg);
    %
    % Create strips of weights, based on the 9-points method with
    % weights 2 for adjacent points and 1 for diagonal points
    wgma = [1 2 1]; 
    nwg = length(wgma);     nma = 2*nwg;
    wgNS = repmat([1 2 1], nrows-1, 1);
    wgWE = repmat([1 2 1]', 1, ncols-1);
    grNS = NaN.*ones(nrows-1,nwg);
    grWE = NaN.*ones(nwg,ncols-1);
    %
    % Calculate time series of spatial gradients from monthly detrended data
    for k = 1:ndata;
        temp = datain(:,:,k);
        % North-South gradients
        diffNS = (temp(1:end-1,:) - temp(2:end,:))/CELLSIZE;
        for j = 1:ncolssg; 
            grNS = wgNS.*diffNS(:,j:j+2);
            sgradNS(:,j,k) = nansum(grNS(1:end-1,:) + grNS(2:end,:),2)/nma;          
        end;
        % West-East gradients 
        diffWE = (temp(:,2:end) - temp(:,1:end-1))/CELLSIZE;
        for i = 1:nrowssg;
            grWE = wgWE.*diffWE(i:i+2,:);
            sgradWE(i,:,k) = nansum(grWE(:,2:end) + grWE(:,1:end-1),1)/nma;
        end;
    end;
    %
    % Initialise work arrays
    tempNS = NaN.*ones(ndata,1);
    tempWE = tempNS;
    tempabs = tempNS;
    % Calculate absolute value (length) of the spatial gradient and
    % complete time series required for the uncertainty analysis
    for m = 1:mdata;
        i = inddata(m,1);
        j = inddata(m,2);
	    tempNS(:) = sgradNS(i,j,:);  
	    tempWE(:) = sgradWE(i,j,:);
	    tempabs = sqrt(tempNS.*tempNS + tempWE.*tempWE); 
        sgradhasdata_abs(:,m) = tempabs;
	    sgradmed_abs(i,j) = median(tempabs);
    end;
    % Clear auxiliary fields
    clear sgradNS sgradWE;
