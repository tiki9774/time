% Time Kalmen Filter


f = [1 10 100 1e3 10e3 100e3]
L = [-50 -70 -113 -128 -135 -140];

sigma = L_F_2_sigma(f,L)

f0 = 10e6;

t = linspace(0,1,1000);

phase = 2*pi*f0*t;
n = sigma .* randn(1,length(t))

phase_n = phase + n;

[p,S] = polyfit(t,phase_n,1);
tmp = polyval(p,t);
figure
plot(t,sin(phase_n),'o')
hold on;
plot(t,sin(tmp))
