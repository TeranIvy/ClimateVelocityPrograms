# ClimateVelocityPrograms

Author: Iva Kavcic

## SaveData

* savedata_serial.m - Reads and saves data in binary *.mat format.

## MannKendall

All the below scripts require function [ktaub.m](http://www.mathworks.com/matlabcentral/fileexchange/authors/23983).

* MKseason_ktaub.m - Mann-Kendall trend analysis on seasonal data.
* MKseason_ktaubMult.m - Mann-Kendall trend analysis on seasonal data
                         with multiple values per season.
* MKyear_ktaub.m - Mann-Kendall trend analysis on year data.

## SpatialGradients

* SpatialGradients_Season.m - Calculates the median spatial gradients.
* SpatGrad_Median.m - Helper function for SpatialGradients_Season.m.
* SpatialGradUncertaintyMatlab_Season.m - Calculates the uncertainties
                                          in median spatial gradients. Before
                                          calculating the spatial gradients the
                                          data are detrended using Sen slopes
                                          from Mann-Kendall trend analysis.
* Detrend_Sen.m - Helper function for SpatialGradUncertaintyMatlab_Season.m.
* SpatGrad_Median_Unc.m - Helper function for SpatialGradUncertaintyMatlab_Season.m.
* SpatGradMatlab_Limits_Unc.m - Helper function for SpatialGradUncertaintyMatlab_Season.m.
