close all
clear all

load res11ev0
ni11ev0angle = ni;

load res11ev45
ni11ev45angle = ni;

load res11ev90
ni11ev90angle = ni;

plot(X/DLe, ni11ev0angle, 'b')
hold
plot(X/DLe, ni11ev45angle, 'r')
plot(X/DLe, ni11ev90angle, 'k')