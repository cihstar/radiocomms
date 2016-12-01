function [ bits ] = DQPSK_demod( s )

% Array to save demodulated bits too
bits = zeros(1,length(s)*2);
bitIndex = 1; %index of this array

for k = 1:length(s)    
    % Save to bit stream
    decision = s(k);
    if(imag(decision) >= 0)
        thisBits(1) = 0;
    else
        thisBits(1) = 1;
    end
    if(real(decision) >= 0)
        thisBits(2) = 0;
    else
        thisBits(2) = 1;
    end
    bits(bitIndex) = thisBits(1);
    bits(bitIndex+1) = thisBits(2);
    bitIndex = bitIndex + 2; 
end
