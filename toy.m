Globals1D;

inB=1.8; outB=11.8;          %x domain [xl,xr]
meshNum = 1000;              
dx = (outB-inB) / meshNum;    %dx: mesh size

x = inB - dx : dx : outB + dx;        % generate the grid point

%Set initial conditions
A = power(x,-2);
B = power(x,-1);
A_exact = A;
B_exact = B;
time = 0;                                                                                                
dt = 0.2*dx;
FinalTime = 10000.0; 
Nsteps = FinalTime/dt;

%Set boundary condition type
DIRICHLET = 1;
FREEZING = 2;
bdry = 1;
                                                                                                         
%list to store the data                                                                              
time_seq = [];
rhs_A_seq = [];
rhs_B_seq = [];

err_A_seq = [];
err_B_seq = [];
                                                                                                         
for tstep=1:Nsteps                                                                                   
    [rhs_A, rhs_B] = RHS_toy(A, B, A_exact, B_exact);
    A = A+rhs_A*dt;                                                                      
    B = B+rhs_B*dt;                                                                      

    % Increment time                                                                                 
    time = time+dt;                                                                                  
    time_seq = [time_seq, time];
    rhs_A_seq = [rhs_A_seq, L2norm(rhs_A)];
    rhs_B_seq = [rhs_B_seq, L2norm(rhs_B)];
    
    err_A_seq = [err_A_seq, L2norm(A-A_exact)];
    err_B_seq = [err_B_seq, L2norm(B-B_exact)];
    if (mod(time, 1) <= dt)
        figure(1); 
        subplot(2,3,1);
        plot(x, A-A_exact); title(['error of A, t = ', num2str(time)]); drawnow;
        subplot(2,3,2);
        plot(x, A); title(['A, t = ', num2str(time)]); drawnow;
        subplot(2,3,3);
        semilogy(time_seq, err_A_seq); title(['L2 norm of error of A with time']); drawnow;
        subplot(2,3,4);
        semilogy(time_seq, err_B_seq); title(['L2 norm of error of B with time']); drawnow;
        subplot(2,3,5);
        semilogy(time_seq, rhs_A_seq); title(['L2 norm of rhs\_A with time']); drawnow; 
        subplot(2,3,6);
        semilogy(time_seq, rhs_B_seq); title(['L2 norm of rhs\_B with time']); drawnow;
    end;
end;  
