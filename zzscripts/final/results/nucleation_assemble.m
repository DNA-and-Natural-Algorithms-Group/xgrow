function [s1,s2,s3,s4,ns1,ns2,ns3,ns4] = nucleation_assemble (filepath)
    %
    % This function assembles the data from seed and no seed files and
    % returns an array of the tao at which each crystal nucleated.
    %
    % each array holds the tao values t which
    %
    s1 = nuc_assemble_one([filepath 'seed/clean/' 'np']);
    s2 = nuc_assemble_one([filepath 'seed/clean/' 'p2']);
    s3 = nuc_assemble_one([filepath 'seed/clean/' 'p3']);
    s4 = nuc_assemble_one([filepath 'seed/clean/' 'p4']);

    ns1 = nuc_assemble_one([filepath 'noseed/clean/' 'np']);
    ns2 = nuc_assemble_one([filepath 'noseed/clean/' 'p2']);
    ns3 = nuc_assemble_one([filepath 'noseed/clean/' 'p3']);
    ns4 = nuc_assemble_one([filepath 'noseed/clean/' 'p4']);

    
% assembles results for one type of proofreading (np, 2p, 3p, 4p)
% only counts flakes larger than 80 tiles. 
% this is necessary because tinybox outputs all flakes, but we
% are only interested in one flake per simulation. 
function [tao] = nuc_assemble_one(filepath)
    % load raw xgrow simulation data
    % raw_data(1,2) is first row, second column.
    raw_data = load(filepath);
    [num_rows,num_col] = size(raw_data);
    
    if (num_rows == 0)
        print('error loading');
    else    
        j=1;
        for (i = 1:num_rows) 
            if (xgrow_row_reader(raw_data(i,:), 'size') >= 80)
                tao(j) = xgrow_row_reader(raw_data(i,:), 'tao');
                j=j+1;
            end
        end
    end
    