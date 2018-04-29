% % % ==================================================================== % % %
% % % 
% % % function Detrend_Sen.m
% % % 
% % % Author: Iva Kavcic
% % % 
% % % Date last modified: 18/05/2017
% % % Date of last comments update: 09/04/2018
% % % Runs with Matlab versions R2010a and newer
% % % 
% % % This function detrends 1 km seasonal data read by the main program
% % % SpatialGradUncertaintyMatlab_Season.m using Sen's slopes calculated
% % % by the Mann-Kendall trend analysis.
% % %
% % % Input variables:
% % %       datain - Input data
% % %       senslope - Sen's slope from Mann-Kendall trend analysis
% % %       timein - Time coordinate from Mann-Kendall trend analysis
% % % Output variables:       
% % %       datadetr - Detrended data
% % %
% % % Sen's slopes were calculated using Jeff Burkey's Matlab script
% % % ktaub.m, which can be retrieved from
% % % http://www.mathworks.com/matlabcentral/fileexchange/authors/23983
% % %  
% % % ==================================================================== % % %
%
function [datadetr] = Detrend_Sen(datain, senslope, timein)  
    %
    % Calculate size of input data
    [nrows, ncols, ndata] = size(datain);
    % Initialise output fields
    datadetr = NaN.*ones(size(datain));
    % Calculate mid time and adjust time coordinate
    midtime = timein(round(length(timein)/2));
    timecoord = timein - midtime;
    % Find median of input data
    meddata = median(datain,3);
    % Detrend using median and adjusted time coordinate
    for k = 1:ndata;
        slope = meddata + senslope.*timecoord(k);
	    datadetr(:,:,k) = datain(:,:,k) - slope;
    end;
    