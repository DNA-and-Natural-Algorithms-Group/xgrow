function nucleation_plot (s1,s2,s3,s4,ns1,ns2,ns3,ns4)
%
% This function takes nucleation data for seed(weave) and no seed and plots
% it in the current open figure.
%

    edges = 1.7:.005:2.1;

    
    
    subplot(411);  
    hold on;  
    xlim([1.7,2.1])
    hist( s1,edges); 
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r')
    hist(ns1,edges);
    
    title('nucleation non-proofreading')
    ylabel('# occurances')
    xlabel('Tao')
    legend({'seed','noseed'})
    
    
    subplot(412);
    hold on;  
    xlim([1.7,2.1])
    hist( s2,edges);  
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r')
    hist(ns2,edges);
    
    title('nucleation 2-redundant-proofreading')
    ylabel('# occurances')
    xlabel('Tao')
    legend({'seed','noseed'})
    
    
    subplot(413);  
    hold on;  
    xlim([1.7,2.1])
    hist( s3,edges);  
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r')
    hist(ns3,edges);
    
    title('nucleation 3-redundant-proofreading')
    ylabel('# occurances')
    xlabel('Tao')
    legend({'seed','noseed'})
    
    
    subplot(414);  
    hold on;  
    xlim([1.7,2.1])
    hist( s4,edges);  
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r')
    hist(ns4,edges);
    
    title('nucleation 4-redundant-proofreading')
    ylabel('# occurances')
    xlabel('Tao')
    legend({'seed','noseed'})
    
    