% % % ==================================================================== % % %
% % % 
% % % function SpatGradMatlab_Limits_Unc.m
% % % 
% % % Author: Iva Kavcic
% % % 
% % % Date last modified: 16/05/2017
% % % Date of last comments update: 09/04/2018
% % % Runs with Matlab versions R2010a and newer
% % % 
% % % This function generates uncertainty limits for the median absolute
% % % spatial gradients of 1 km seasonal data. The uncertainty limits
% % % are generated using Matlab functions bootsrp and bootci on the
% % % time series of detrended absolute spatial gradients. The time
% % % series are processed in chunks for faster performance.
% % % The function returns uncertainty limits and spread to the 
% % % program SpatialGradUncertaintyMatlab_Season.m.
% % %
% % % Input variables:
% % %       sgraddata - Reshaped time series of absolute spatial
% % %                   gradients for the uncertainty analysis, 
% % %                   dimension [ndata, mdata]
% % %       nrowssg - Number of rows of 2D median spatial gradient field
% % %       ncolssg - Number of columns of 2D median spatial gradient
% % %                 field
% % %       inddata - Row and column indices of the complete spatial 
% % %                 gradient time series, dimension [mdata,2]
% % % Output variables:   
% % %       sgraddata_limits - Uncertainty limits for medians of absolute
% % %                          spatial gradients, based on 5 and 95
% % %                          percentiles of the boostrapped time series,
% % %                          dimension [nrowssg, ncolssg, 2]
% % %       sgraddata_stdev - Spread (1 stdev) for medians of absolute
% % %                         spatial gradients, dimension
% % %                         [nrowssg, ncolssg, 2]
% % %
% % % Spatial gradients on detrended data were calculated in the 
% % % function SpatGrad_Median_Unc.m using the 9-points method from
% % % Burrows et al. 2011.
% % %  
% % % ==================================================================== % % %
%
function [sgraddata_limits, sgraddata_stdev] = ...
          SpatGrad_Limits_Unc(sgraddata, nrowssg, ncolssg, inddata)
    %
    % Calculate size of input data and allocate output data
    [ndata, mdata] = size(sgraddata);
    sgraddata_limits = NaN.*ones(nrowssg, ncolssg, 2);
    sgraddata_stdev = NaN.*ones(nrowssg, ncolssg);
    %
    % Calculate size of blocks to process
    chunksize = 1000
    indc(:,1) = [1:chunksize:mdata];
    indc(:,2) = [chunksize:chunksize:mdata mdata];
    nchunks = size(indc,1);
    %
    % Generate 1000 boostrapped time series for each complete time limits
    % of spatial gradients and calculate uncertainty limits and spread
    nreps = 1000;
    for k = 1:nchunks;
        % Set start and end indices of a chunk
        k1 = indc(k,1); k2 = indc(k,2);
        display([num2str(k1) ', ' num2str(k2)]);
        offset = (k-1)*chunksize;
        nelem = k2 - k1 + 1;
        % Store the chunk in a temporary array
        dZ = sgraddata(:,k1:k2);
        % Calculate spread of the median spatial gradients
        stdtemp(1:nelem) = std(bootstrp(nreps,@median,dZ));
        % Calculate uncertainty limits of the median spatial gradients
        climtemp(1:2,1:nelem) = bootci(nreps,{@median,dZ},'alpha',0.05);
        % Store results for each row and column index in the chunk
        for m = k1:k2;      
            i = inddata(m,1);
            j = inddata(m,2);
            sgraddata_stdev(i, j) = stdtemp(m-offset);
            sgraddata_limits(i, j, 1) = climtemp(1,m-offset);
            sgraddata_limits(i, j, 2) = climtemp(2,m-offset);
        end;      
    end;
   