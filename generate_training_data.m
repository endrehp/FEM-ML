%Generate training data

Omegas = [10,15,17,20,50];

UNL_total = [];
F_total =[];

counter = 1;
for i=1:5

    tf = 1 + randi([1,3]);
    
   OMEGA = Omegas(counter);
   
   Nonlinear_Cantilever_Vibration_varying_input
   
   UNL_total = [UNL_total, UNL]; 
   F_total = [F_total, FNLin];
   
   if counter == length(Omegas)
       counter = 1;
   end
    
   counter
end
