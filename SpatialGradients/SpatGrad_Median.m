% % % ==================================================================== % % %
% % % 
% % % function SpatGrad_Median.m
% % % 
% % % Author: Iva Kavcic
% % % 
% % % Date last modified: 19/05/2017
% % % Date of last comments update: 08/04/2018
% % % Runs with Matlab versions R2010a and newer
% % % 
% % % This function uses 1 km seasonal data read by the main program
% % % SpatialGradients_Season.m and to calculate spatial gradients of data.
% % % It returns medians of spatial gradients in NS, WE and abs directions,
% % % as well as their angles.
% % %
% % % Input variables:
% % %       datain - Input data
% % %       CELLSIZE - Size of spatial cell
% % % Output variables:       
% % %       sgradmed_NS - Spatial gradients in North-South direction,
% % %                     medians over years
% % %       sgradmed_WE - Spatial gradients in West-East direction,
% % %                     medians over years   
% % %       sgradmed_abs - Absolute values of spatial gradients,
% % %                      medians over years 
% % %       sgradmed_angle - Angles of median spatial gradients     
% % %       inddata - Indices of rows and columns for complete time series,
% % %                 dimension [mdata,2]
% % %
% % % Spatial gradients are calculated by the 9-points method described in
% % % Burrows et al. 2011, The Pace of Shifting Climate in Marine and
% % % Terrestrial Ecosystems (Science 334:652-655),
% % % http://science.sciencemag.org/content/334/6056/652
% % %  
% % % ==================================================================== % % %
%
function [sgradmed_NS, sgradmed_WE, sgradmed_abs, sgradmed_angle, ...
          inddata] = SpatGrad_Median(datain, CELLSIZE)  
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
    % Find row and column indices with complete time series
    [inddata(:,1), inddata(:,2)] = find(~isnan(hasalldata)); 
    mdata = size(inddata,1);
    %
	% Initialise output and auxiliary fields
    sgradNS = NaN.*ones(nrowssg,ncolssg,ndata);
    sgradWE = sgradNS;
    sgradhasdata_abs = NaN.*ones(ndata,mdata);
    sgradmed_NS = NaN.*ones(nrowssg,ncolssg);
    sgradmed_WE = sgradmed_NS;
    sgradmed_abs = sgradmed_NS;
    sgradmed_angle = sgradmed_NS;
    %
    % Create strips of weights, based on the 9-points method with
    % weights 2 for adjacent points and 1 for diagonal points
    wgma = [1 2 1]; 
    nwg = length(wgma);    nma = 2*nwg;
    wgNS = repmat([1 2 1], nrows-1, 1);
    wgWE = repmat([1 2 1]', 1, ncols-1);
    grNS = NaN.*ones(nrows-1,nwg);
    grWE = NaN.*ones(nwg,ncols-1);
    %
    % Calculate time series of spatial gradients from monthly data
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
    % Calculate absolute value (length) of the spatial gradient, relevant
    % medians and angles for the complete time series
    for m = 1:mdata;
        i = inddata(m,1);
        j = inddata(m,2);
	    tempNS(:) = sgradNS(i,j,:);  
	    tempWE(:) = sgradWE(i,j,:);
	    tempabs = sqrt(tempNS.*tempNS + tempWE.*tempWE);
	    sgradmed_NS(i,j) = median(tempNS);
	    sgradmed_WE(i,j) = median(tempWE);
	    sgradmed_abs(i,j) = median(tempabs);
	    % Spatial gradient angles in degrees, range [-180,180]
	    sgradmed_angle(i,j) = atan2(sgradmed_NS(i,j),sgradmed_WE(i,j)).*180./pi;
    end;
    % Clear auxiliary fields
    clear sgradNS sgradWE;
