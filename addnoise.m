function y = addnoise(x,SNR)

%   Adds noise to signal x with predefined SNR

noise = randn(size(x));
Noise_Power = norm(noise,2);
Signal_Power = norm(x-mean(x),2);
% Initial_SNR = 20*(log10(Signal_Power/Noise_Power))

K = (Signal_Power/Noise_Power)^2*10^(-SNR/10);  % Scale factor

y = x + sqrt(K)*noise; 

end

