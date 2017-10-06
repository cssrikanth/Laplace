clc
clear all
close all
data = 'runtime.dat';
m = dlmread(data);
sec=m(:,2);

figure(1)
bar(sec);
xlabel ("Processors");
ylabel ("Average time(sec)");
title ("Scalability for a 1000 x 1000 2D grid");
saveas (1, "scalability.eps")