
snrs = 0:0.5:10;
bers_dqpsk_orig = zeros(1,length(snrs));
bers_dqpsk = zeros(1,length(snrs));
bers_dqpsk_b = zeros(1,length(snrs));
bers_theo = zeros(1,length(snrs));
ber_index = 1;
ber_total = 0;
N = 3;
for SNR = snrs
    ber_total = 0;
    for n = 1:N
        ber_total = ber_total + DigiTrans_with_IF(SNR,0,1,1);
    end
    bers_dqpsk(ber_index) = ber_total/N;
    
    ber_total = 0;
    for n = 1:N
        ber_total = ber_total + DigiTrans_with_IF(SNR,0,1,0);
    end
    bers_dqpsk_b(ber_index) = ber_total/N;
    
    EbNo_lin = 10^(SNR/10);
    bers_theo(ber_index) = 0.5 * erfc(sqrt(2*EbNo_lin)/sqrt(2));
    ber_index = ber_index + 1;
end

figure(2);

semilogy(snrs,bers_theo,'r',snrs,bers_dqpsk_b,'bx', snrs,bers_dqpsk,'gx','LineWidth',2);
legend('Theoretical QSPK', 'Simulated baseband DQPSK', 'Simulated IF DQPSK');
title('BER v SNR');
xlabel('SNR / dB');
ylabel('BER');
grid on;