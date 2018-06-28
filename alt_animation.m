close all

%d_max = max(max(UNL_total));

d=100;
x = linspace(0,1, length(predictions(1,1:d:end)));


for i=1:d
    
    cla
    clf
   
    p_d = predictions(1, i:d:end); %dependent predictions
    %p_i = preds(i, :); %independent predictions
    y_t = y_test(1, i:d:end); %test data
    %ul = UL_total(1:2:end, end-length(preds(:,1))+i)/d_max; %linear solution
    
    
    hold on
    plot(x, p_d, '-o', 'DisplayName', 'dependent')
    %plot(x, p_i, '-', 'DisplayName', 'independent')
    plot(x, y_t, 'DisplayName', 'test')
    %plot(x, ul, '-x')
    %plot((unl-ul).^2)
    %disp(max((unl-ul).^2))
    axis([0 1.2 -1 1])
    %legend('show')
    drawnow
    pause(0.1)
    

end