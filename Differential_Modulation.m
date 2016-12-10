function s_hat = Differential_Modulation(s)
    delay = 1;
    s_hat = zeros(1,length(s));
    for i = 1:length(s)
        s_hat(i) = delay * s(i);
        delay = s_hat(i);
    end
end

