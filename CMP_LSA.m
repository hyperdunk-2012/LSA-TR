%% LSA-TR
I = 2;

SymbolTR = zeros(K,num);
SymbolTRAIN = zeros(K,num);

tic

for m = 1:num
	T = zeros(K,1);
	S = SymbolTX(1:K,m);
    s = TX(1:K,m);
    f = zeros(K,1);
    s0 = zeros(K,1);

    for n = 1:I

        Signal_Power = abs(s); 
        Mean_Power   = mean(Signal_Power);
        M = 2*Mean_Power;
        
        for k = 1:K
            if abs(s(k,1)) >= M
                s0(k,1) = M*exp(1i*phase(s(k,1)));
            else
                s0(k,1) = s(k,1);
            end
        end
        
        f = s0-s;
        F = fft(f)/sqrt(K);

        for k = 1:L
            T(TR_pos(k),1) = F(TR_pos(k),1);
        end
       
        t = sqrt(K)*ifft(T);
        
        P1 = 0;
        P2 = 0;
        for k = 1:K
            if f(k,1) ~= 0
                P1 = P1+abs(f(k,1))*abs(t(k,1));
                P2 = P2+abs(t(k,1))*abs(t(k,1));
            end
        end
        P = P1/P2;
        
        s = s+P*t;

    end
    
    TR(1:K,m) = s;
    
end

toc
tim=toc;

Signal_Power = abs(TR.^2);
Peak_Power   = max(Signal_Power,[],1);
Mean_Power   = mean(Signal_Power,1);
PAPR_TR3 = 10*log10(Peak_Power./Mean_Power);

[cdf3, PAPR3] = ecdf(PAPR_TR3);
ecdf3 = 1-cdf3;

figure(3)
semilogy(PAPR1,ecdf1,'-b',PAPR3,ecdf3,'-r')
ylim([0.001,1]);
legend('Original','LSA-TR')
title('COMPARE')
xlabel('PAPR0 [dB]');
ylabel('CCDF (Pr[PAPR>PAPR0])');
grid on
save(filename);