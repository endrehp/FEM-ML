close all

x = linspace(0,1, length(preds(1,:)));
for i=1:length(preds(:,1))
    
    cla
    clf
   
    p_d = predictions(i, :); %dependent predictions
    p_i = preds(i, :); %independent predictions
    y_t = y_test(i,:); %test data
    
    
    hold on
    plot(x, p_d, 'o', 'DisplayName', 'dependent')
    plot(x, p_i, '-', 'DisplayName', 'independent')
    plot(x, y_t, 'DisplayName', 'test')
    %plot((unl-ul).^2)
    %disp(max((unl-ul).^2))
    axis([0 1.2 -1 1])
    %legend('show')
    drawnow
    pause(0.01)
    

end