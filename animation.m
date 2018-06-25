%drawsoltion

size_y = size(y);
n_timesteps = size_y(1);

for i=1:size(y)
    cla;
    clf;
    dz_lin =  zeros(10,1);
    dz_nl =  zeros(10,1);
    dz_lin(2:end) = y_lin(i,1:2:17);
    dz_nl(2:end) = y_nl(i,1:2:17);
    
    hold on
    plot(x,dz_lin, 'o')
    plot(x, dz_nl)
    axis([0 1 0 40]);
    drawnow
end