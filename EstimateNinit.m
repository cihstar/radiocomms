function Ninit = EstimateNinit(signal, N, num_symbols, mu0, beta, d0)
mu = mu0 * (N/4);

signal_length = length(signal);
delay = d0 * ones(1,round(signal_length/N)); 
index=1; 
early = 0;
late = 0;

%moving_average_max = 100;
%moving_average = zeros(1,moving_average_max);
%moving_average_index = 1;

for k = (N+1): N : (N*num_symbols)
    lastDelay = round(delay(index)*N);    
    late  = (sqrt(1-beta)^2 * signal(k + lastDelay + 1)) - (beta * late);
    early = (sqrt(1-beta)^2 * signal(k + lastDelay - 1)) - (beta * early);
    delay(index+1) = delay(index) + mu * (abs(late)^2 - abs(early)^2); 
    
    %{
    moving_average(moving_average_index) = delay(index+1)*N;    
    moving_average_index = moving_average_index + 1;
    if(moving_average_index == moving_average_max + 1)
        moving_average_index = 1;
    end
    
    if(mod(index, moving_average_max) == 0)
        if(max(moving_average) - min(moving_average) < 0.4)
            Ninit = round(mean(moving_average)) + 1;
            break
        end
    end
    %}
    index = index+1;
end
%index
plot((delay(1:index)*N)+1,'k');
Ninit = mode(round((delay(1:index)*N)))+1;
end

      
            
        
        

