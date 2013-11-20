function [fn, sn] = powspec( deltaT, yn )

nT = size(yn,1);
maxT = nT*deltaT;

nyquistFrequencyT = 1/(2*deltaT);	% nyquist frequency
fourierFrequencyT = 1/(nT*deltaT);	% fourier frequency
fn = ([0:ceil(nT/2)-1 -floor(nT/2):-1]*fourierFrequencyT)';

ybar = fft(yn)*deltaT/sqrt(maxT);

sn = ybar .* conj(ybar);

fn = fftshift(fn);
sn = fftshift(sn,1);