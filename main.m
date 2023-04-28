% Time Kalmen Filter

clear all;
close all;


f = [.1 1 10 100 1e3 10e3 100e3];
L = [-30 -50 -70 -113 -128 -135 -140];


sigma = L_F_2_sigma(f,L)   %radians


%%
f0 = 1e6;
fs = f0*2.5;
t = 0:1/fs:.1;

phase = 2*pi*f0*t;
n = sigma .* randn(1,length(t));

phase_n = phase + n;

[p,S] = polyfit(t,phase_n,1);
tmp = polyval(p,t);
% figure
% plot(t,sin(phase_n),'o')
% hold on;
% plot(t,sin(tmp))

v = sin(phase_n);
N = length(v);
xdft = fft(v);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(v):fs/2;
figure
plot(freq,pow2db(psdx))
grid on
title("Periodogram Using FFT")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")




