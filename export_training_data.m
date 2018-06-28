close all

%Export training data as csv

%FM = [F_total; UNL_total];
filename = 'training_data.csv';
csvwrite(filename,FM')