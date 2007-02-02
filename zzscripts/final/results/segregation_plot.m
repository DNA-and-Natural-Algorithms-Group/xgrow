function segregation_plot()
    label = ['0000 ' '1111 ' '1010 ' '0100 ' '1000 ' '1100 ' '1001 ' '0110 ' '1110 ' '1101 '];
    expected = [0.06 0.06 0.13 0.13 0.13 0.13 0.06 0.06 0.13 0.13];
    noweave = [0.47 0.16 0.05 0.06 0.08 0.07 0.04 0.02 0.02 0.03];
    weave = [0.26 0.02 0.23 0.17 0.11 0.07 0.06 0.04 0.02 0.02];
    
    hold on
    bar([0 1 2 3 4 5 6 7 8 9],[expected(:), noweave(:), weave(:)], 'grouped')
    xlabel(['pattern: ' label])
    ylabel('proportion observed')
end