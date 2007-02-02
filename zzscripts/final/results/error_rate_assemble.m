function [np_err, p2_err, p3_err, p4_err] = error_rates_assemble(filepath)
%
% This function assembles the data from np 2p 3p and 4p data files and 
% returns the error rates vs tao of each.
%
% *_err is a 2 dimensional array
% -- dim1 : tao (gmc/gse) (unknown # levels)
% -- dim2 : error_rate (1), error_bar (2)

    np_err = assemble_one([filepath 'pn'],1,4);
    p2_err = assemble_one([filepath 'p2'],2,2);
    p3_err = assemble_one([filepath 'p3'],3,1);
    p4_err = assemble_one([filepath 'p4'],4,1);
end

% does the actual work described above.
% --loads data file
%   --loops through each line : taos
%   --loops through each line : ave error rate for each tao
%   --loops through each line : std dev of error rate
function [errors] = assemble_one(filepath, n, x)
    % load raw xgrow simulation data
    raw_data = load(filepath);
    [num_rows,num_col] = size(raw_data);
    
    [taos, num_taos] = tally_taos(raw_data, num_rows);
    taos = sort(taos);
    [errs,num_samples] = tally_error(raw_data, num_rows, ...
        taos, num_taos,  n, x);
    [stddevs] = tally_stddev(raw_data, num_rows, taos, num_taos, ...
           n, x, errs, num_samples);
 
    errors = cleanup(taos,errs,stddevs,num_taos,num_samples);
end

% return an aray of all tao values in data
% and the total number of tao values
function [taos, num_taos] = tally_taos(raw_data, num_rows)
    % count current number of tao values
    num_taos = 0;
    for (i = 1:num_rows)
        curr_tao = xgrow_row_reader(raw_data(i,:), 'tao');
        
        if (num_taos > 0 && tao_row(curr_tao, taos(:), num_taos))
            % do nothing, this row is alreayd in the array of taos
        else
            % add current value to the array of taos
            num_taos = num_taos+1;
            taos(num_taos) = curr_tao; 
        end
    end
end

% determine the average error rate for each tao
function [error_rate, num_samples] = tally_error(raw_data, num_rows, ...
                                                taos, num_taos, n, x)
    num_samples = zeros(1,num_taos);
    error_rate = zeros(1,num_taos);
    for (i = 1:num_rows)
        curr_tao = xgrow_row_reader(raw_data(i,:), 'tao');
        
        row=tao_row(curr_tao, taos(:), num_taos);
        % determine this row's error rate
        err = xgrow_row_reader(raw_data(i,:), 'mismatches') / n;
        att = xgrow_row_reader(raw_data(i,:), 'rows') * x;
        temp = err/att;
        error_rate(row) = error_rate(row) + temp;

        % keep running count of the number of items in each tao bin
        num_samples(row) = num_samples(row) + 1;
    end
    
    % divide each bin by # items in bin to get the average
    error_rate = error_rate ./ num_samples;
end

% returns twice the standard deviation of the sample for each tao
function [stddev] = tally_stddev(raw_data, num_rows, ...
                                 taos, num_taos, n, x, ...
                                 av_error, num_samples)
    stddev =  zeros(1,num_taos);
    for (i = 1:num_rows)
        curr_tao = xgrow_row_reader(raw_data(i,:), 'tao');
        
        row = tao_row(curr_tao, taos(:), num_taos);
        
        % determine this row's error rate
        err = xgrow_row_reader(raw_data(i,:), 'mismatches') / n;
        att = xgrow_row_reader(raw_data(i,:), 'rows') * x;
        temp = err/att;
        
        % determine this row's deviation
        dev = (temp - av_error(row))^2;
        stddev(row) = stddev(row) + dev;
    end
    
    % do stddev touchups
    stddev = stddev ./ num_samples;
    stddev = sqrt(stddev);
    stddev = 2 * stddev;
end

% determine which row the current tao is in, returns false
% if tao is not in array yet.
function [result] = tao_row(tao, all_taos, num_taos)
   result = false;
   for(i = 1:num_taos)
       if (tao == all_taos(i))
           result = i;
       end
   end
end

function [result] = cleanup(taos,errs,stds,num_taos,num_samples)
    j = 0;
    for (i = 1:num_taos)
        if (num_samples(i) > 50)
            j = j + 1;
            result(1,j) = taos(i);
            result(2,j) = errs(i);
            result(3,j) = stds(i);
        end
    end
end