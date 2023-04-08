% Time Kalmen Filter
fc = 10 %Hz
t = linspace(0,1,10000)
x1 = sin(2*pi*fc*t);
x2 = cos(2*pi*fc*t);

x3 = x1 + x2

figure;
plot(x1)
hold on;
plot(x2)
hold on;
plot(x3)