function exp_plot()
    label = ['np(no weave)' 'np control' 'np' '2p' '3p' '4p'];
    err = [0.033997804 0.007887911 0.009652676 0.003579638 0.001565753 0.000988443];
    err_bar = [0.010533846 0.001372893 0.001207845 0.000860332 0.001182671 0.000883653];
    
    hold on
    bar(err);
    errorbar(err, err_bar, '.');
    xlabel('experment')
    ylabel('error rate')
end
