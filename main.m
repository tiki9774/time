% Time Kalmen Filter

clear all;
close all;
clc;

f = [.1 1 10 100 1e3 10e3 100e3];
L = [-30 -50 -70 -113 -128 -135 -140];


sigma = L_F_2_sigma(f,L)   %radians


%%
f0 = 1e6;
fs = f0*5;
t = 0:1/fs:1;

phase = 2*pi*f0*t;
n = sigma .* randn(1,length(t));
phase_n = zeros(1,length(t));
dt = 1/fs;
for ii = 2:length(t)
    phase_n(ii) = phase_n(ii-1) + 2*pi*f0*dt+(n(ii));
end
[p,S] = polyfit(t,phase_n,1);
tmp = polyval(p,t);
% figure
% plot(t,sin(phase_n),'o')
% hold on;
% plot(t,sin(tmp))

figure
plot(phase_n-phase)

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

[M,I] = max(psdx);
psdx = psdx/ M;
psdx = psdx(I:end);
freq = freq(I-1) - freq(I:end);
figure;
plot(10*log10(freq),10*log10(psdx))
grid on;

% %%
% clear all;
% close all;
% 
% f0 = 1;
% fs = f0*5;
% t = 0:1/fs:1000000;
% phi0 = 0;
% n = randn(1,length(t));
% dt = 1/fs;
% 
% phase = zeros(1,length(t));
% phase(1) = phi0;
%dt = 1/fs
% for ii = 2:length(t)
%     phase(ii) = phase(ii-1) + 2*pi*f0*dt + n(ii);
% end
% 
% [p,S] = polyfit(t,phase,1);
% tmp = polyval(p,t);
% figure
% plot(t,phase,'o')
% hold on;
% plot(t,tmp)
% 
% delta = phase - tmp;
% figure
% plot(delta)


