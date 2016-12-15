function [ s ] = QPSK_mod( bits )

B = (2:2:length(bits));
s = zeros(1,length(bits)/2);
for k = B
    if bits(k-1) == 1 && bits(k) == 1
        s(k/2) = 1 * exp(1j*-(3/4)*pi);        
    elseif bits(k-1) == 0 && bits(k) == 1
        s(k/2) = 1 * exp(1j*(3/4)*pi);   
    elseif bits(k-1) == 0 && bits(k) == 0
        s(k/2) = 1 * exp(1j*(1/4)*pi);       
    elseif bits(k-1) == 1 && bits(k) == 0
        s(k/2) = 1 * exp(1j*-(1/4)*pi);
    end
end