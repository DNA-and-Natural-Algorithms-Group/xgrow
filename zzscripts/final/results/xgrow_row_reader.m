function [result] = xgrow_row_reader(row, attribute)
%
% this file provides functions that obtain data from an xgrow result file
% row format : (this is the default for an xgrow file)
%  1    2   3   4      5        6        7        8       9      10
% gmc, gse, k, time, size, mismatches, events, perimeter, dG, dG_bonds 
%
% enter a row of data, and a string attribute (the ones shown above)
% and it will return the value of that attribute for this row.
%
% ie xgrow_row_reader(row, 'gmc') will return the gmc of that row.
%
% additional attributes: 
%      tao = gmc/gse
%      rows = size / 6                  
    if (strcmpi(attribute, 'gmc'))        loc = 1; end   
    if (strcmpi(attribute, 'gse'))        loc = 2; end
    if (strcmpi(attribute, 'size'))       loc = 5; end
    if (strcmpi(attribute, 'mismatches')) loc = 6; end
        
    if (strcmpi(attribute, 'tao')) 
        % tao = gmc / gse
        result = row(1)/row(2);
        return;
    end

    if (strcmpi(attribute, 'rows'))
        result = row(5) / 6;
        return;
    end
        
    result = row(loc);
    