Globals1D;
GlobalsGR;

inB=1.8; outB=10;          %x domain [xl,xr]
meshNum = 1000;              
dx = (outB-inB) / meshNum;    %dx: mesh size
dt = 0.2*dx;

% Evaluate the initial conditions
x = inB - dx : dx : outB + dx;        % generate the grid point

% compute time step size                                                                                 
time = 0;                                                                                                
FinalTime = 10000.0; 
Nsteps = FinalTime/dt;

printf("N=%d\n", 0);
printf("mesh numbers=%d\n", meshNum);

% Set initial conditions
%init_min;
%init_iso;
%init_PG_out;
%init_PG_in;
%init_ks_in;
init_ks_out;

%Set bdry condition
DIRICHLET = 1;
FREEZING = 2;
bdry_type = 1; 

%Evolve S or psi
USE_S = 0;
USE_psi = 0;

%lists to store results
time_seq = [];
rhs_g11_seq = [];
rhs_Pi11_seq = [];
rhs_Phi11_seq = [];

err_g11_seq = [];
err_Pi11_seq = [];
err_Phi11_seq = [];

C0_seq = [];
C1_seq = [];
Cr00_seq = [];
Cr01_seq = [];
Cr11_seq = [];
AH_seq = [];
                                                                                                         
for tstep=1:Nsteps                                                                                   
    compute_RHS;
    g00 = g00+dt*rhs_g00;                
    g01 = g01+dt*rhs_g01;                
    g11 = g11+dt*rhs_g11;                
    Pi00 = Pi00+dt*rhs_Pi00;             
    Pi01 = Pi01+dt*rhs_Pi01;             
    Pi11 = Pi11+dt*rhs_Pi11;             
    Phi00 = Phi00+dt*rhs_Phi00;          
    Phi01 = Phi01+dt*rhs_Phi01;          
    Phi11 = Phi11+dt*rhs_Phi11;          
    if (USE_S)
        S = S+dt*rhs_S;                  
        Pi_S = Pi_S+dt*rhs_Pi_S;         
        Phi_S = Phi_S+dt*rhs_Phi_S;      
    end;
    if (USE_psi)
        psi = psi+dt*rhs_psi;            
        Pi_psi = Pi_psi+dt*rhs_Pi_psi;   
        Phi_psi = Phi_psi+dt*rhs_Phi_psi;
    end;

    if (bdry_type == FREEZING)
        apply_freezing;
    end;

    % Increment time                                                                                 
    time = time+dt;                                                                                  

    time_seq = [time_seq, time];
    rhs_g11_seq = [rhs_g11_seq, L2norm(rhs_g11)];
    rhs_Pi11_seq = [rhs_Pi11_seq, L2norm(rhs_Pi11)];
    rhs_Phi11_seq = [rhs_Phi11_seq, L2norm(rhs_Phi11)];

    err_g11_seq = [err_g11_seq, L2norm(g11-g11_exact)];
    err_Pi11_seq = [err_Pi11_seq, L2norm(Pi11-Pi11_exact)];
    err_Phi11_seq = [err_Phi11_seq, L2norm(Phi11-Phi11_exact)];

    C0_seq = [C0_seq, L2norm(C0)];
    C1_seq = [C1_seq, L2norm(C1)];
    Cr00_seq = [Cr00_seq, L2norm(Cr00)];
    Cr01_seq = [Cr01_seq, L2norm(Cr01)];
    Cr11_seq = [Cr11_seq, L2norm(Cr11)];
    for (i=1:size(AH_indicator)(2)-1)
        if (AH_indicator(i) <= 0.0 && AH_indicator(i+1) > 0.0)
            AH_seq = [AH_seq, (x(i+1)*AH_indicator(i)-x(i)*AH_indicator(i+1))/(AH_indicator(i)-AH_indicator(i+1)) - 2.0];
            break;
        end
    end
    
    if (mod(time, 5) < dt)
        figure(1);
        subplot(2,4,1);
        plot(x, g11-g11_exact, 'r'); title(['Error of g11, Pi11, Phi11, t = ', num2str(time)]); hold on;
        plot(x, Pi11-Pi11_exact, 'b'); hold on;
        plot(x, Phi11-Phi11_exact, 'y'); hold off;
        subplot(2,4,2);
        semilogy(time_seq, err_g11_seq, 'r'); title(['L2 norm of error with time']); hold on;
        semilogy(time_seq, err_Pi11_seq, 'b'); hold on; 
        semilogy(time_seq, err_Phi11_seq, 'y'); hold off;
        subplot(2,4,3);
        plot(x, C0, 'r'); title(['C0 and C1, t = ', num2str(time)]); hold on; plot(x, C1, 'b'); hold off; subplot(2,4,4);
        semilogy(time_seq, C0_seq, 'r'); title(['L2 norm of C0 and C1 with time']); hold on;
        semilogy(time_seq, C1_seq, 'b'); hold off;
        subplot(2,4,5);
        plot(x, g11, 'r'); title(['g11, Pi11, Phi11, t = ', num2str(time)]); hold on;
        plot(x, Pi11, 'b'); hold on;
        plot(x, Phi11, 'y'); hold off;
        %plot(x, Cr00, 'r'); title(['Cr, t= ', num2str(time)]); hold on;
        %plot(x, Cr01, 'b'); hold on;
        %plot(x, Cr11, 'y'); hold off;
        subplot(2,4,6);
        plot(x, Cr00, 'r'); title(['Cr, t= ', num2str(time)]); hold on;
        plot(x, Cr01, 'b'); hold on;
        plot(x, Cr11, 'y'); hold off;

        %semilogy(time_seq, Cr00_seq, 'r'); title(['L2 norm of Cr with time']); hold on;
        %semilogy(time_seq, Cr01_seq, 'b'); hold on;
        %semilogy(time_seq, Cr11_seq, 'y'); hold off;

        %plot(x, rhs_g11, 'r'); title(['rhs, t= ', num2str(time)]); hold on;
        %plot(x, rhs_Pi11, 'b'); hold on;
        %plot(x, rhs_Phi11, 'y'); hold off;

        subplot(2,4,7);
        semilogy(time_seq, rhs_g11_seq, 'r'); title(['L2 norm of rhs with time']); hold on;
        semilogy(time_seq, rhs_Pi11_seq, 'b'); hold on;
        semilogy(time_seq, rhs_Phi11_seq, 'y'); hold off;

        subplot(2,4,8);
        %plot(x, k1); title(['dS\_1/dg11']); drawnow;
        plot(x, AH_indicator); title(['AH\_indicator, t = ', num2str(time)]); drawnow;
        %plot(time_seq, AH_seq, 'r.'); title(['AH position drift with time']); drawnow;

    end;
    if (time >= 1000 && time - dt <= 1000)
        printf("L2 norm of error of g11 when t=20:  %.2e\n", err_g11_seq(end));
        printf("L2 norm of error of Pi11 when t=20:  %.2e\n", err_Pi11_seq(end));
        printf("L2 norm of error of Phi11 when t=20:  %.2e\n", err_Phi11_seq(end));
        break;
    end
end;  
