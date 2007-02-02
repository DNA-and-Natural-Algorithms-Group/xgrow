%
% This function takes error rate data for np, 2p 3p, 4p and plots it in the
% current open figure. converts error rate data to log scale before
% plotting
%
% see error_rate_assemble for rquired format of error_rates array
%
function error_rate_plot (np_err, p2_err, p3_err, p4_err) 
  % subplot(121) 
%    semilogy( np_err(1,:),np_err(2,:), 's:y', ...
%              p2_err(1,:),p2_err(2,:), '+:g', ...
%              p3_err(1,:),p3_err(2,:), 'o:b', ...
%              p4_err(1,:),p4_err(2,:), 'x:m')
   hold on
   err_bar(np_err, '-y.');
   err_bar(p2_err, '-g.');
   err_bar(p3_err, '-b.');
   err_bar(p4_err, '-m.');
   
   legend({'np','2p','3p','4p'})
   title('simulation error rates')
   xlabel('Tao')
   ylabel('error rate')
   set(gca,'XGrid','on')
   set(gca,'YGrid','on')
         
  % subplot(122)
  % hold on
  % lin(np_err, p2_err, p3_err, p4_err)
end

function err_bar(data, color)
    a = log(data(2,:));
    up = log(data(2,:) + data(3,:)) - a;
    low = a - log(data(2,:) - data(3,:));
    errorbar(data(1,:), a,  up, color)
end

function lin(np_err, p2_err, p3_err, p4_err)
    set(gca,'XGrid','on')
    set(gca,'YGrid','on')
    errorbar(np_err(1,:), np_err(2,:), np_err(3,:), 'y')
    errorbar(p2_err(1,:), p2_err(2,:), p2_err(3,:), 'g')
    errorbar(p3_err(1,:), p3_err(2,:), p3_err(3,:), 'b')
    errorbar(p4_err(1,:), p4_err(2,:), p4_err(3,:), 'm')
end