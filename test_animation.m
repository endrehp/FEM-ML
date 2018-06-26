close all

x = linspace(0,1, length(preds(1,:)));
d_max = max(max(UNL_total));
for i=1:length(preds(:,1))
    
    cla
    clf
   
    p_d = predictions(i, :); %dependent predictions
    p_i = preds(i, :); %independent predictions
    y_t = y_test(i,:); %test data
    %ul = UL_total(1:2:end, end-length(preds(:,1))+i)/d_max; %linear solution
    
    
    hold on
    plot(x, p_d, '-o', 'DisplayName', 'dependent')
    %plot(x, p_i, '-', 'DisplayName', 'independent')
    plot(x, y_t, 'DisplayName', 'test')
    %plot(x, ul, '-x')
    %plot((unl-ul).^2)
    %disp(max((unl-ul).^2))
    axis([0 1.2 -0.4 0.4])
    %legend('show')
    drawnow
    pause(0.01)
    

end