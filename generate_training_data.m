clear all
close all
%Generate training data

Omegas = [16.5,17,18,20,22,17.68];

UNL_total = [];
UL_total = [];
F_total =[];

counter_h = 1;
for i=1:5

    tf = 1 + randi([1,3]);
    
   OMEGA = Omegas(counter_h);
   
   Nonlinear_Cantilever_Vibration_varying_input
   
   UNL_total = [UNL_total, UNL]; 
   UL_total = [UL_total, ULin];
   F_total = [F_total, FNLin];
   
   if counter_h == length(Omegas)
       counter_h = 1;
   end
    
   counter_h = counter_h + 1;
   counter_h
end
