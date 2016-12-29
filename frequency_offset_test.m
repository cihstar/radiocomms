%%script to plot BER curves for varying frequency offset
frequency = [-99 0 99];
SNR = 0:13;
figure()
hold on
for f = frequency
    BER = generate_BER_plot(SNR, 'DigiTrans_baseband', [0 1 f]);
    semilogy(SNR, BER)
end
legend(strread(num2str(frequency),'%s'))
ax = gca
ax.YLabel.String = 'BER'
ax.XLabel.String = 'SNR'
ax.YScale = 'log'

    