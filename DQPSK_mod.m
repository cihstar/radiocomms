function [ s ] = DQPSK_mod( bits )

B = (2:2:length(bits));
delta_phase_shifts = zeros(1,length(bits)/2);
phase_index = 1;
for k = B
    if bits(k-1) == 1 && bits(k) == 1
        delta_phase_shifts(phase_index) = -(3/4)*pi;        
    elseif bits(k-1) == 0 && bits(k) == 1
        delta_phase_shifts(phase_index) = (3/4)*pi;        
    elseif bits(k-1) == 0 && bits(k) == 0
        delta_phase_shifts(phase_index) = (1/4)*pi;        
    elseif bits(k-1) == 1 && bits(k) == 0
        delta_phase_shifts(phase_index) = -(1/4)*pi;
    end
    phase_index = phase_index + 1;
end

s = zeros(1,length(delta_phase_shifts));
s(1) = 1;
for k = 2:length(s)
    s(k) = s(k-1) * exp(1i*delta_phase_shifts(k-1));
end