function   [Ap, As, Ai] = OPOsol(Ap,As,Ai,kappa_p,kappa_s,kappa_i,dk,Z2,dzLN)
        
        %dk = ks + ki - kp;
        %k1 = f(Xn,Yn)
        kp_1 = (1j*kappa_p)*(As.*Ai)*exp(1j*dk*Z2);        %k1 for pump
        ks_1 = (1j*kappa_s)*(conj(Ai).*Ap)*exp(-1j*dk*Z2);   %k1 for signal
        ki_1 = (1j*kappa_i)*(conj(As).*Ap)*exp(-1j*dk*Z2);   %k1 for idler
             
        %A_half2 = A + k1*dz/2;
        Ap_half2 = Ap + kp_1*dzLN/2; 
        As_half2 = As + ks_1*dzLN/2; 
        Ai_half2 = Ai + ki_1*dzLN/2; 
        %k2 = f(Xn+h/2,Yn+k1*h/2) 
        kp_2 = (1j*kappa_p)*(As_half2.*Ai_half2)*exp(1j*dk*(Z2+dzLN/2)); %k2 for pump
        ks_2 = (1j*kappa_s)*(conj(Ai_half2).*Ap_half2)*exp(-1j*dk*(Z2+dzLN/2));%k2 for signal
        ki_2 = (1j*kappa_i)*(conj(As_half2).*Ap_half2)*exp(-1j*dk*(Z2+dzLN/2));     %k2 for idler
        
        %A_half3 = A + k2*dz/2;
        Ap_half3 = Ap + kp_2*dzLN/2; 
        As_half3 = As + ks_2*dzLN/2; 
        Ai_half3 = Ai + ki_2*dzLN/2; 
        %k3 = f(Xn+h/2,Yn+k2*h/2)  
        kp_3 = (1j*kappa_p)*(As_half3.*Ai_half3)*exp(1j*dk*(Z2+dzLN/2)); %k3 for pump
        ks_3 = (1j*kappa_s)*(conj(Ai_half3).*Ap_half3)*exp(-1j*dk*(Z2+dzLN/2));%k3 for signal
        ki_3 = (1j*kappa_i)*(conj(As_half3).*Ap_half3)*exp(-1j*dk*(Z2+dzLN/2));      %k3 for idler
      
        %A_full = A + k3*dz;
        Ap_full = Ap + kp_3*dzLN; 
        As_full = As + ks_3*dzLN; 
        Ai_full = Ai + ki_3*dzLN; 
        %k4 = f(Xn+h,Yn+k3*h)
        kp_4 = (1j*kappa_p)*(As_full.*Ai_full)*exp(1j*dk*(Z2+dzLN));      %k4 for pump
        ks_4 = (1j*kappa_s)*(conj(Ai_full).*Ap_full)*exp(-1j*dk*(Z2+dzLN));   %k4 for signal
        ki_4 = (1j*kappa_i)*(conj(As_full).*Ap_full)*exp(-1j*dk*(Z2+dzLN));         %k4 for idler
        
        %A = A + dz*(k1 + 2*k2 + 2*k3 + k4)/6;
        Ap = Ap + dzLN*(kp_1 + 2*kp_2 + 2*kp_3 + kp_4)/6;
        As = As + dzLN*(ks_1 + 2*ks_2 + 2*ks_3 + ks_4)/6;
        Ai = Ai + dzLN*(ki_1 + 2*ki_2 + 2*ki_3 + ki_4)/6;
        
end