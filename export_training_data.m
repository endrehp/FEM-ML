close all

%Export training data as csv

FM = [FNLin; UNL];
filename = 'training_data.csv';
csvwrite(filename,FM')