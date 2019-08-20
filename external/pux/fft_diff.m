function [zp,zpp] = fft_diff(z,idx)
%Computes the first and second derivatives of a number of drop boundaries
%using the FFT.
zp = zeros(size(z));
zpp = zeros(size(z));
for j = 1:size(idx,1)
    N = idx(j,2)-idx(j,1)+1;
    c = fft(z(idx(j,1):idx(j,2)));
    k = [(0:N/2) (-N/2+1:-1)]';
    tmp = c.*k;
    zp(idx(j,1):idx(j,2)) = ifft(1i*tmp);
    zpp(idx(j,1):idx(j,2)) = ifft(-tmp.*k);
end

