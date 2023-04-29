function sigma_phasenoise = L_F_2_sigma(f,L)
%inputs
% L: L(f) in dBc/Hz 
% f: frequency in Hz
%outputs
% sigma_phasenoise: std of phasenoise of L(f) in rad

%

% model phasenoise conversion as sum of piecewise functions
% each region in peicewise function is trapazoid sum
if length(f) ~= length(L)
    error('Input vector lengths dont match')
end

L_f = 10 .^ ( L ./ 10);
phasenoise_contributions = zeros(length(f)-1,1);
for ii= 1:(length(f)-1)
    phasenoise_contributions(ii) = (f(ii+1) - f(ii)) * (L_f(ii+1) + L_f(ii));
end



sigma_phasenoise = sqrt(sum(phasenoise_contributions));
end