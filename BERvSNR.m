snrs = 0:0.25:10;
bers = zeros(1,length(snrs));
bers_theo = zeros(1,length(snrs));
ber_index = 1;
ber_total = 0;
N = 20;
for SNR = snrs
    ber_total = 0;
    for n = 1:N
        ber_total = ber_total + DigiTrans_baseband(SNR,0,0);
    end
    bers(ber_index) = ber_total/N;
    bers_theo(ber_index) = qfunc(sqrt(2*SNR));
    ber_index = ber_index + 1;
end

figure(1);
semilogy(snrs,bers,'b',snrs,bers_theo,'r');
title('BER v SNR');
xlabel('SNR / dB');
ylabel('BER');