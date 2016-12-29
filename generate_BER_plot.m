function [BER_p] = generate_BER_plot(SNR_arr, fun, args)
%% fun = string name of function, e.g. 'DigiTrans_baseband' with the first 
%           input argument being SNR, other input arguments defined in args

% minBER = ths function will keep reducing SNR until min BER is reached
% args = other inputs to fun, not including SNR

i = 1;
SNR_arr;
for SNR = SNR_arr
    in = [{fun} {SNR}, num2cell(args)];
    BER_p(i) = feval(in{:});
    i = i+1;
end

end


