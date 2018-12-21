clc;
clear all;
close all;
syms t tau % Short-cut for constructing symbolic objects.
x = heaviside(t-3)-heaviside(t-5); % input signal to the LTI system
h = exp(-3*t)*heaviside(t); % system’s impulse response
y=int(subs(x,tau)*subs(h,t-tau),tau,0,t)
disp('The output of the convolution integral is');
disp(y);
subplot(3,1,1);
ezplot(x,[0,10]);

title('Input signal for LTI system');
xlabel('t')
ylabel('x(t)')
grid;
subplot(3,1,2);
ezplot(h);

title('LTI System impulse response');
xlabel('t')
ylabel('h(t)')
grid;
ti = 0;
tf = 10;
N = 1000;
dt = (tf-ti)/N;
t = ti:dt:tf;
x = heaviside(t-3)-heaviside(t-5);
x=x*dt;
if t>=0
    h = exp(-3*t);
else
    h=0;
end
y1=conv(x,h);
subplot(3,1,3);
plot(t, y1 (1: length (t)));
title('Output of LTI system');
xlabel('t')
ylabel('y(t)')
grid;