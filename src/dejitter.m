function f = dejitter(im_in, M, par_alpha)
%DEJITTER Summary of this function goes here
%   Detailed explanation goes here
[r,c] = size(im_in);
N = M + 1;

ft = padarray(im_in(:,1),N,'both');
f = padarray(im_in,[0 N],'both');

gamma = im_in(:,(N+1):c-N);

p = zeros(1,r);

p(1) = N+1;
p(2) = p(1);

phi1 = padarray(gamma(1,:)',N,'both');
phi2 = phi1;

h = zeros(c,2*N+1);

gad = zeros(1,2*N+1);

for i = 2:r
    for k = 1:(2*N+1)
        h(:,k) = [zeros(1,k-1) gamma(i,:) zeros(1,2*N - k + 1)];
        m = max([k,p(i),p(i-1)]);
        n = min([k,p(i),p(i-1)]) + c - (2*N + 1);
        n
        
        gad(k) = 1/(n - m + 1) * sum(abs(h(m:n,k) - 2*phi1(m:n) +phi2(m:n)) .^ par_alpha);
        
    end
    
    [min_v,min_i] = min(gad);
    p(i+1) = min_i;
    
    phi2 = phi1;
    phi1 = h(:,min_i);
    
    f(:,i) = [zeros(1,min_i - 1) im_in(i,:) zeros(1,2*N - min_i + 1)];
    
end






end

