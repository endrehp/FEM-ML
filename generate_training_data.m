%clear all
close all
%Generate training data

%Omegas = [16.5,17,18,20,22,17.68];

X_train_total = [];
UL_total = [];
Y_train_total =[];
OMEGA = 6;
f_constant = 0.0001;
counter_h = 1;
for i=1:10

    %tf = 1 + randi([1,3]);
    tf = 3;
    
   %OMEGA = Omegas(counter_h);
   
   Nonlinear_Cantilever_Vibration_varying_input
   
   FM = FNLin';
   DM = UNL';
   
   dim = size(DM);
   
   F = FM(1:10:end,1:2:end);
   F = F/0.0001;
    D = DM(1:10:end,1:2:end);
    D = D/0.03;
    
    [n_timesteps, n_nodes] = size(D);
    d = 200; %number of steps in "model memory"
    t = d;
    
    FL = zeros(n_timesteps+d, length(F(1,:)));
    FL(d+1:end,:) = F;
    DL = zeros(n_timesteps+d, length(D(1,:)));
    DL(d+1:end,:) = D;
   
    
    X_train = zeros(n_timesteps, (d+1)*n_nodes);
    Y_train = zeros(n_timesteps, n_nodes);
    
    for i=1:(n_timesteps-1)
    
    bulk = 1;
    for j=1:(n_nodes)-1
        X_train(i, bulk:bulk + d) = FL(t-d+1:t+1,j);    
        bulk = bulk + d+1;
    
    end
    Y_train(i, :) = DL(t, :);
    t = t + 1;
    
    end
    
    X_train_total = [X_train_total; X_train];
    Y_train_total = [Y_train_total; Y_train];
   
    
   f_constant = 0.0001*(-1+2*rand());
   counter_h = counter_h + 1;
   counter_h
   
end

%Export to csv

x_filename = 'x_train.csv';
y_filename = 'y_train.csv';

csvwrite(x_filename,X_train_total)
csvwrite(y_filename,Y_train_total)
