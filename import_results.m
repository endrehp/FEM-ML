close all

%Import results

predictions = csvread('predictions.csv');
preds = csvread('preds.csv');
y_test = csvread('y_test.csv');
y_train = csvread('y_train.csv');
