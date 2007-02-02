function zzsim_results()
    
    % SURF 2006 ZigZag redundancy proofreading
    % script created Jan 2007
    %
    % Simulation results for zigzag proofreading
    % -- error rates
    % -- nucleation A (gmc matches error rate gmc =13)
    % -- nucleation B (gmc realistic =???)

    %include('error_rate_plot');
    %include('nucleation_plot');
    %include('nucleation_assemble');
    %include('error_rate_assemble');

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% error rates %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

    % assembling data from files
      [np_err, p2_err, p3_err, p4_err] = error_rate_assemble('../errorrates/clean/');

    % plotting data
      figure(1)
      error_rate_plot(np_err, p2_err, p3_err, p4_err);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% nucleation A %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
    
    % assmebling data from files
      [s1,s2,s3,s4,  ns1,ns2,ns3,ns4] = nucleation_assemble('../nucA/');

    % plotting data
      figure(2)
      nucleation_plot(s1,s2,s3,s4,  ns1,ns2,ns3,ns4);

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% nuclation B %%%%%
%%%%%%%%%%%%%%%%%%%%%%%
    
    % assembling data from files
      %[s1,s2,s3,s4,  ns1,ns2,ns3,ns4] = nucleation_assemble('../nucB/');

    % plotting data
      %figure(3)
      %nucleation_plot(s1,s2,s3,s4,  ns1,ns2,ns3,ns4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% experimental results %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % pattern proportion
    figure(4)
    segregation_plot();
    
    % error rates
    figure(5)
    exp_plot();