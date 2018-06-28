# FEM-ML

This is a summer research project where the goal is to use an artificial neural network to predict the response of a nonlinear dynamic problem.
Nonlinear dynamic FEA can be time consuming, this can be a problem if we want to do simulations in real time. Therefore we will try to train 
an ANN on a large number of load cases to make it able to predict the response of a variety of loads in a small fraction of the time it 
would take to solve the fully nonlinear FEM problem again. This approach will hopefully allow us to run a real time simulation of the dynamic
system.

At the moment the training data is created by a nonlinear dynamic beam model in Matlab.
