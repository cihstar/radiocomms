function [s_hat] = DQPSK_diff_detection(s, N)
    s_hat = zeros(1,length(s));
    s_hat(1:N) = s(1:N);
    for x = N : N : length(s)-1
        for y = 1:N
            s_hat(x+y) = s(x+y) * conj( s( (x-N)+y ) );
        end
    end
end