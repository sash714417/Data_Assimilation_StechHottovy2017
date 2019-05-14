function [A_k,rlv,qeqs] = SWEModelAk(modes,k,tau1,tau2,b1,b2,tildeQ1,tildeQ2,taur,taul,tauv)

    mcheck = 5*modes+3;

    rlv = zeros(3*(modes+1),mcheck);
    num_q0 = 3*(modes+1);
    
    % Equations
    % r0
    rlv(1,1) = -1i*k-1/taur;
    rlv(1,num_q0+1) = -1/sqrt(2)*1/tau1;
    rlv(1,num_q0+2) = -1/sqrt(2)*1/tau2;
    
    % r1
    rlv(2,2) = -1i*k-1/taur;
    rlv(2,3) = 1;
    rlv(2,num_q0+3) =-1/sqrt(2)*1/tau1;
    rlv(2,num_q0+4) =-1/sqrt(2)*1/tau2;
    % v0
    rlv(3,2) = -1;
    rlv(3,3) = -1/tauv; 

    for j=3:(modes+2)
        
        % Define equation number
        num_eq = 3*(j-2)+1;

        % rn
        rlv(num_eq,num_eq) = -1i*k-1/taur;
        rlv(num_eq,num_eq+2) = sqrt(j-1);
        if j-1 < modes
            rlv(num_eq,num_q0+2*(j-1)+1) = -1/sqrt(2)*1/tau1;
            rlv(num_eq,num_q0+2*(j-1)+2) = -1/sqrt(2)*1/tau2;
        end
        
        % ln-2
        rlv(num_eq+1,num_eq+1) = 1i*k-1/taul;
        rlv(num_eq+1,num_eq+2) = sqrt(j-2);
        rlv(num_eq+1,num_q0+2*(j-3)+1) = 1/sqrt(2)*1/tau1;
        rlv(num_eq+1,num_q0+2*(j-3)+2) = 1/sqrt(2)*1/tau2;
        
        % vn-1
        rlv(num_eq+2,num_eq)=-sqrt(j-1);
        rlv(num_eq+2,num_eq+1) = -sqrt(j-2);
        rlv(num_eq+2,num_eq+2) = -1/tauv;
               
    end
    
    % Build q equations
    qeqs = zeros(2*modes,mcheck);
    
    %q0 equation
    
    %r0 and l0 contributions
    qeqs(1,1) = -tildeQ1*1i*k/sqrt(2);
    qeqs(2,1) = -tildeQ2*1i*k/sqrt(2);
    
    qeqs(1,5) = -tildeQ1*1i*k/sqrt(2);
    qeqs(2,5) = -tildeQ2*1i*k/sqrt(2);
    
    % v contributions
    qeqs(1,6) = -tildeQ1/sqrt(2);
    qeqs(2,6) = -tildeQ2/sqrt(2);
    
    %q contributions
    qeqs(1,num_q0+1) = -1/tau1-b1*k^2-b1*.5;
    qeqs(2,num_q0+2) = -1/tau2-b2*k^2-b2*.5;
    qeqs(1,num_q0+5) = sqrt(2)*b1/2;
    qeqs(2,num_q0+6) = sqrt(2)*b2/2;
    
    %q1 equation
    
    %r1 and l1 contributions
    qeqs(3,2) = -tildeQ1*1i*k/sqrt(2);
    qeqs(4,2) = -tildeQ2*1i*k/sqrt(2);
    
    qeqs(3,8) = -tildeQ1*1i*k/sqrt(2);
    qeqs(4,8) = -tildeQ2*1i*k/sqrt(2);
   
    %v contributions
    qeqs(3,3) = tildeQ1/sqrt(2);
    qeqs(4,3) = tildeQ2/sqrt(2);
    
    qeqs(3,9) = -tildeQ1;
    qeqs(4,9) = -tildeQ2;

    
    %q contributions
    qeqs(3,num_q0+3) = -1/tau1-b1*k^2-b1*3/2;
    qeqs(4,num_q0+4) = -1/tau2-b2*k^2-b2*3/2;
    
    if modes>3
        qeqs(3,num_q0+7) = sqrt(6)*b1/2;
        qeqs(4,num_q0+8) = sqrt(6)*b2/2;
    end
    
    % Cycle through modes
    for j = 3:modes           

       % rj contributions
       qeqs(2*j-1,3*(j-2)+1) = -tildeQ1*1i*k/sqrt(2);
       qeqs(2*j,3*(j-2)+1) = -tildeQ2*1i*k/sqrt(2);
       
       % vj contributions
       qeqs(2*j-1,3*(j-1)) = tildeQ1*sqrt(j-1)/sqrt(2);
       qeqs(2*j-1,3*(j-1)+6) = -tildeQ1*sqrt(j)/sqrt(2);
       
       qeqs(2*j,3*(j-1)) = tildeQ2*sqrt((j-1))/sqrt(2);
       qeqs(2*j,3*(j-1)+6) = -tildeQ2*sqrt(j)/sqrt(2);
       
       % lj contributions
       qeqs(2*j-1,3*j+2) = -tildeQ1*1i*k/sqrt(2);
       qeqs(2*j,3*j+2) = -tildeQ2*1i*k/sqrt(2);
       
       
       %q contributions
       qeqs(2*j-1,num_q0+2*(j-1)+1) = -1/tau1-b1*k^2-b1*(2*j-1)/2;
      % [tau1, b1, k, (2*j-1)/2]
       qeqs(2*j-1,num_q0+2*(j-3)+1) = b1*sqrt((j-2)*(j-1))/2;
       
       qeqs(2*j,num_q0+2*(j-1)+2) = -1/tau2-b2*k^2-b2*(2*j-1)/2;
       qeqs(2*j,num_q0+2*(j-3)+2) = b2*sqrt((j-2)*(j-1))/2; 
       
       if j+2<modes+1
           qeqs(2*j-1,num_q0+2*(j+1)+1) = b1*sqrt((j)*(j+1))/2;
           qeqs(2*j,num_q0+2*(j+1)+2) = b2*sqrt((j)*(j+1))/2;
       end
        
    end
    
    A_k = [rlv;qeqs];
    
end