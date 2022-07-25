function apply_freezing
Globals1D;
GlobalsGR;

%inner boundary
if (inB>2)
rhs_g00_in = rhs_g00(1);
rhs_g01_in = rhs_g01(1);
rhs_g11_in = rhs_g11(1);
 
rhs_Pi00_in = rhs_Pi00(1);
rhs_Pi01_in = rhs_Pi01(1);
rhs_Pi11_in = rhs_Pi11(1);
                
rhs_Phi00_in = rhs_Phi00(1);
rhs_Phi01_in = rhs_Phi01(1);
rhs_Phi11_in = rhs_Phi11(1);

b = sqrt(1/g11(1));
rhs_g00(1) = 0.;
rhs_g01(1) = 0.;
rhs_g11(1) = 0.;

rhs_Pi00(1) = -paragamma2/2*rhs_g00_in + 1/2*rhs_Pi00_in + b/2*rhs_Phi00_in;
rhs_Pi01(1) = -paragamma2/2*rhs_g01_in + 1/2*rhs_Pi01_in + b/2*rhs_Phi01_in;
rhs_Pi11(1) = -paragamma2/2*rhs_g11_in + 1/2*rhs_Pi11_in + b/2*rhs_Phi11_in;
 
rhs_Phi00(1) = -paragamma2/2/b*rhs_g00_in + 1/2/b*rhs_Pi00_in + 1/2*rhs_Phi00_in  ;
rhs_Phi01(1) = -paragamma2/2/b*rhs_g01_in + 1/2/b*rhs_Pi01_in + 1/2*rhs_Phi01_in  ;
rhs_Phi11(1) = -paragamma2/2/b*rhs_g11_in + 1/2/b*rhs_Pi11_in + 1/2*rhs_Phi11_in  ;
end;

%out bdry
rhs_g00_out = rhs_g00(end);
rhs_g01_out = rhs_g01(end);
rhs_g11_out = rhs_g11(end);
 
rhs_Pi00_out = rhs_Pi00(end);
rhs_Pi01_out = rhs_Pi01(end);
rhs_Pi11_out = rhs_Pi11(end);
                
rhs_Phi00_out = rhs_Phi00(end);
rhs_Phi01_out = rhs_Phi01(end);
rhs_Phi11_out = rhs_Phi11(end);

%rhs_S_out = rhs_S(end);
%rhs_Pi_S_out = rhs_Pi_S(end);
%rhs_Phi_S_out = rhs_Phi_S(end);

b = sqrt(1/g11(end));
%rhs_g00(end) = 0.;
%rhs_g01(end) = 0.;
%rhs_g11(end) = 0.;

rhs_Pi00(end) = paragamma2/2*rhs_g00_out + 1/2*rhs_Pi00_out ...
        + b/2*rhs_Phi00_out;% - 0.25*(Phi00(end)-Phi00_exact(end))*power(b,3)*rhs_g11_out;
rhs_Pi01(end) = paragamma2/2*rhs_g01_out + 1/2*rhs_Pi01_out ...
        + b/2*rhs_Phi01_out;% - 0.25*(Phi01(end)-Phi01_exact(end))*power(b,3)*rhs_g11_out;
rhs_Pi11(end) = paragamma2/2*rhs_g11_out + 1/2*rhs_Pi11_out ...
        + b/2*rhs_Phi11_out;% - 0.25*(Phi11(end)-Phi11_exact(end))*power(b,3)*rhs_g11_out;
rhs_Phi00(end) = -paragamma2/2/b*rhs_g00_out ...
        + 1/2/b*rhs_Pi00_out + 1/2*rhs_Phi00_out;% + 0.25*(Phi00(end)-Phi00_exact(end))*power(b,2)*rhs_g11_out;
rhs_Phi01(end) = -paragamma2/2/b*rhs_g01_out ...
        + 1/2/b*rhs_Pi01_out + 1/2*rhs_Phi01_out;% + 0.25*(Phi01(end)-Phi01_exact(end))*power(b,2)*rhs_g11_out;
rhs_Phi11(end) = -paragamma2/2/b*rhs_g11_out ...
        + 1/2/b*rhs_Pi11_out + 1/2*rhs_Phi11_out;% + 0.25*(Phi11(end)-Phi11_exact(end))*power(b,2)*rhs_g11_out;

%rhs_S(end) = 0.;
%rhs_Pi_S(end) = -paragamma2/2*rhs_S_out + 1/2*rhs_Pi_S_out ...
%        + b/2*rhs_Phi_S_out ;
%rhs_Pi_psi(end) = paragamma2/2*rhs_psi_out + 1/2*rhs_Pi_psi_out ...
%        + b/2*rhs_Phi_psi_out ;
%
%rhs_Phi_S(end) = -paragamma2/2/b*rhs_S_out ...
%        + 1/2/b*rhs_Pi_S_out + 1/2*rhs_Phi_S_out ;
%rhs_Phi_psi(end) = -paragamma2/2/b*rhs_psi_out ...
%        + 1/2/b*rhs_Pi_psi_out + 1/2*rhs_Phi_psi_out ;

%rhs_Pi00(end) = 0.;
%rhs_Pi01(end) = 0.; 
%rhs_Pi11(end) = 0.; 
%rhs_Phi00(end) = 0.;
%rhs_Phi01(end) = 0.;
%rhs_Phi11(end) = 0.;

%rhs_S = 0.*x;
%rhs_Pi_S = 0.*x;
%rhs_Phi_S = 0.*x;
%rhs_psi = 0.*x;
%rhs_Pi_psi = 0.*x;
%rhs_Phi_psi = 0.*x;
