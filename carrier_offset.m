function [ out ] = carrier_offset( in, offset, N )
%in is an array containing the values outputed by a transceiver
%offset: carrier frequency offset as a fraction of the symbol frequency
%N: number of elements per symbol

len = length(in);

for iter = 1:len
    out(iter) = in(iter) * exp(i * (offset* (1/N) * iter) * 2 * pi);
end

end

