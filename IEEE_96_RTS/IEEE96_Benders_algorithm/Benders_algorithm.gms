$ontext
Setting Reserve Requirements to Approximate The Stochastic Dispatch Efficiency
Copyright (C) 2018 Vladimir Dvorkin

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
$offtext

$eolcom //
option iterlim=999999999;          // avoid limit on iterations
option reslim=100000;               // timelimit for solver in sec.
option optcr=0;                     // gap tolerance
option sysout = off;
option solprint=OFF;                // include solution print in .lst file
option limrow=0;                    // limit number of rows in .lst file
option limcol=0;                    // limit number of columns in .lst file
option decimals=1;                  // defines decimals in output file
option mip=cplex;

$Include data_import.gms

*** Multicut Bender's algorithm
************ Initialize the Benders algorithm ************
Sets
         nu       'set of interations'       /nu1*nu8/;
Alias(nu,rho)
;
*** Algorithm-related parameters of the master problem
Parameters
Save_C_RT(s,rho)
Save_alpha(s,rho)
* save subproblem dual values
Save_theta_R_U(i,s,rho)
Save_theta_R_D(i,s,rho)
Save_theta_P_W(k,s,rho)
Save_theta_delta_DA(n,s,rho)
* save day_ahead dispatch
Save_P(i,rho)
Save_DD_U(rho,z,t)
Save_DD_D(rho,z,t)
Save_LB(rho,t)
Save_UB(rho,t)

Save_R_U(i,rho)
Save_R_D(i,rho)
Save_P_W(k,rho)
Save_delta_DA(n,rho)
* Time specific parameters
Parameter
D(j)
W_exp(k)

* iteration counter
Scalar iter;
************ Master problem **************
* Primal variables
Positive variables p(i),p_w(k),R_U(i),R_D(i),DD_RU(z),DD_RD(z);
Free variables Z_M, delta_d(n);
delta_d.FX('n1') =0;
*Dual variables
Free variables alpha(s),lambda_RU,lambda_RD,lambda_D(n);
*reserve clearing
Positive variables rho_U(i),rho_D(i),eta_U(i),eta_D(i),gamma(i),
*energy clearing
mu_U(i),mu_L(i),phi_U(k),phi_L(k),mu_delta(n,m);
$ontext
SOS1 variables
*reserve clearing
u_U(i,y),u_D(i,y),u_U_plus(i,y),u_D_plus(i,y),u_UD(i,y),
*energy clearing
u_under(i,y),u_over(i,y),u_w_under(k,y),u_w_over(k,y),u_delta(n,m,y);
$offtext
Binary variables
*reserve clearing
u_U(i),u_D(i),u_U_plus(i),u_D_plus(i),u_UD(i),
*energy clearing
u_under(i),u_over(i),u_w_under(k),u_w_over(k),u_delta(n,m);
Equations
OBJ_M,CUT_zero,CUT
*Stationarity conditions
STAT1(i),STAT2(i),STAT3(i),STAT4(k),STAT5(n),
*Primary feasibility
*reserve
URB,DRB,URLd(i),DRLd(i),UDB(i),
*energy
PB_da(n),MaPC_c(i),MiPC_c(i),MPC_w(k),MTC(n,m),
*Complimentarity slakness
*reserve
CS11(i),CS12(i),CS21(i),CS22(i),CS31(i),CS32(i),CS41(i),CS42(i),CS101(i),CS102(i),
*energy
CS51(i),CS52(i),CS61(i),CS62(i),CS71(k),CS72(k),CS81(k),CS82(k),CS91(n,m),CS92(n,m)
;
OBJ_M..
         Z_M =e= sum(i, C(i)*p(i) + C_up(i)*R_U(i) + C_dw(i)*R_D(i)) + sum(s, pi(s) * alpha(s));
* Cut constraints
CUT_zero(s)..
         alpha(s) =G= -10e7;
CUT(rho,s)$(iter>1 and ord(rho)<iter)..
         alpha(s) =G= Save_C_RT(s,rho)
         + sum((i), Save_theta_R_U(i,s,rho) * [R_U(i)  - Save_R_U(i,rho)])
         + sum((i), Save_theta_R_D(i,s,rho) * [R_D(i)  - Save_R_D(i,rho)])
         + sum((k), Save_theta_P_W(k,s,rho) * [P_W(k)  - Save_P_W(k,rho)])
         + sum((n), Save_theta_delta_DA(n,s,rho) * [delta_d(n) - Save_delta_DA(n,rho)])
         ;
STAT1(i)..
         C_up(i) - lambda_RU + gamma(i) + rho_U(i) - eta_U(i) =e= 0;
STAT2(i)..
         C_dw(i) - lambda_RD + gamma(i) + rho_D(i) - eta_D(i) =e= 0;
STAT3(i)..
         c(i) + mu_U(i) - mu_L(i) - sum(n$MG(i,n), lambda_D(n)) =E= 0;
STAT4(k)..
         phi_U(k) - phi_L(k) - sum(n$MW(k,n), lambda_D(n)) =E= 0;
STAT5(n)..
         sum(m$(MN(n,m)=1), (mu_delta(n,m)+lambda_D(n))/X(n,m) - (mu_delta(m,n)+lambda_D(m))/X(m,n)) =E= 0;
*Primary feasibility
*reserve
URB(z)..
         sum(i$MZ(i,z), R_U(i)) =e= DD_RU(z);
DRB(z)..
         sum(i$MZ(i,z), R_D(i)) =e= DD_RD(z);
URLd(i)..
         R_U(i) =l= R_up_max(i);
DRLd(i)..
         R_D(i) =l= R_dw_max(i);
UDB(i)..
         R_U(i) + R_D(i) =L= P_max(i);
*energy
PB_da(n)..
         sum(i$MG(i,n), p(i)) + sum(k$MW(k,n), p_w(k))  - sum(j$ML(j,n), D(j))
         =e=
         sum(m$(MN(n,m)=1),(delta_d(n)-delta_d(m))/X(n,m));
MaPC_c(i)..
         p(i) =l= P_max(i) - R_U(i);
MiPC_c(i)..
         p(i) =g= R_D(i);
MPC_w(k)..
         p_w(k) =l= W_exp(k);
MTC(n,m)$MN(n,m)..
         (1/X(n,m))*(delta_d(n)-delta_d(m)) =l= T_cap(n,m);
*Complimentarity slakness
*reserve
CS11(i)..
         rho_U(i) =l= MM * u_U(i);
CS12(i)..
         R_up_max(i) - R_U(i) =l= R_up_max(i) * (1-u_U(i));
CS21(i)..
         rho_D(i) =l= MM * u_D(i);
CS22(i)..
         R_dw_max(i) - R_D(i) =l= R_dw_max(i) * (1-u_D(i));
CS31(i)..
         eta_U(i) =l= MM * u_U_plus(i);
CS32(i)..
         R_U(i) =l= R_up_max(i) * (1-u_U_plus(i));
CS41(i)..
         eta_D(i) =l= MM * u_D_plus(i);
CS42(i)..
         R_D(i) =l= R_dw_max(i) * (1-u_D_plus(i));
CS101(i)..
         gamma(i) =L= MM * u_UD(i);
CS102(i)..
         P_max(i) - R_U(i)- R_D(i) =L= P_max(i) *(1- u_UD(i));
*energy
CS51(i)..
         mu_U(i) =L= MM * u_over(i);
CS52(i)..
         P_max(i)- R_U(i) - p(i) =L= P_max(i) * (1-u_over(i));
CS61(i)..
         mu_L(i) =L= MM * u_under(i);
CS62(i)..
         - R_D(i) + p(i) =L= P_max(i)*(1-u_under(i));
CS71(k)..
         phi_U(k) =L= MM * u_w_over(k);
CS72(k)..
         W_exp(k) - p_w(k) =L= W_exp(k) * (1-u_w_over(k));
CS81(k)..
         phi_L(k) =L= MM * u_w_under(k);
CS82(k)..
         p_w(k) =L= W_exp(k) * (1-u_w_under(k));
CS91(n,m)$(MN(n,m)=1)..
         mu_delta(n,m) =L= MM * u_delta(n,m);
CS92(n,m)$(MN(n,m)=1)..
         T_cap(n,m) - (1/X(n,m))*(delta_d(n)-delta_d(m)) =L= 2*T_cap(n,m)* (1 - u_delta(n,m));
Model MASTER /all/;

************ Subproblems problem **************
*** Algorithm-related parameters of the subproblems
Parameters P_a_s(k),FIX_R_U(i),FIX_R_D(i),FIX_P_W(k),FIX_delta_DA(n);

Free variables Z_S, delta_o(n);
Positive variables rr_u(i),rr_d(i),w_sp(k),l(j);
delta_o.FX('n1') = 0;
Equations OBJ_S,PB_o(n),URDL(i),DRDL(i),MTC_o(n,m),WSL(k),SHL(j),SENS1(i),SENS2(i),SENS3(k),SENS4(n);
OBJ_S..
         Z_S =E= sum(i, C_pl(i)*rr_u(i) - C_mn(i)*rr_d(i)) + sum(j, VOLL*l(j));
PB_o(n)..
         sum(i$MG(i,n), rr_u(i)-rr_d(i)) + sum(j$ML(j,n), l(j))
         + sum(k$MW(k,n), P_a_s(k)) - sum(k$MW(k,n), p_w(k))
         - sum(k$MW(k,n), w_sp(k)) =e=
         sum(m$(MN(n,m)=1),(1/X(n,m))*(delta_o(n)-delta_d(n)-delta_o(m)+delta_d(m)));
URDL(i)..
         rr_u(i) =l= R_U(i);
DRDL(i)..
         rr_d(i) =l= R_D(i);
MTC_o(n,m)$MN(n,m) ..
         (1/X(n,m))*(delta_o(n)-delta_o(m)) =l= T_cap(n,m);
WSL(k)..
         w_sp(k) =l= P_a_s(k);
SHL(j)..
         l(j) =l= D(j);
SENS1(i)..
         R_U(i) =E= FIX_R_U(i);
SENS2(i)..
         R_D(i) =E= FIX_R_D(i);
SENS3(k)..
         p_w(k) =E= FIX_P_W(k);
SENS4(n)..
         delta_d(n) =E= FIX_delta_DA(n);
Model SUBPROBLEM /OBJ_S,PB_o,URDL,DRDL,MTC_o,WSL,SHL,SENS1,SENS2,SENS3,SENS4/;

************ Benders' Algorithm ************
* Initialize parameters
Parameters
UB(nu)
LB(nu)
;

Scalar conv /10e3/;
Scalar epsilon /0.01/;

Save_C_RT(s,rho) = 0;
Save_alpha(s,rho) = 0;

Save_theta_R_U(i,s,rho) = 0;
Save_theta_R_D(i,s,rho) = 0;
Save_theta_P_W(k,s,rho) = 0;
Save_theta_delta_DA(n,s,rho) = 0;

Save_R_U(i,rho) = 0;
Save_R_D(i,rho) = 0;
Save_P_W(k,rho) = 0;
Save_delta_DA(n,rho) = 0;

Parameter P_a_t(k,s);
Parameter pen;
pen = 0.25;

Loop(t,
*IF(ord(t)>7 and ord(t)<24,
* Data
D(j)=D_data(j,t);
P_a_t(k,s)=P_b(s,t)*wind_cap(k)*pen;
W_exp(k)=sum(s, pi(s)*P_a_t(k,s));
Loop(nu$(conv ge epsilon),
         iter = ord(nu);

         Solve MASTER using mip minimizing Z_M;

         Save_R_U(i,nu) = R_U.l(i);
         Save_R_D(i,nu) = R_D.l(i);
         Save_P_W(k,nu) = P_W.l(k);
         Save_delta_DA(n,nu) = delta_d.l(n);

         Save_P(i,nu) = P.l(i);
         Save_DD_U(nu,z,t) = DD_RU.l(z);
         Save_DD_D(nu,z,t) = DD_RD.l(z);
         Save_alpha(s,nu) = alpha.l(s);

         LB(nu) = Z_M.l;

         Loop(s,
                  P_a_s(k) = P_a_t(k,s);
                  FIX_R_U(i) = Save_R_U(i,nu);
                  FIX_R_D(i) = Save_R_D(i,nu);
                  FIX_P_W(k) = Save_P_W(k,nu);
                  FIX_delta_DA(n) = Save_delta_DA(n,nu);

                  Solve SUBPROBLEM using mip minimizing Z_S;

                  Save_theta_R_U(i,s,nu) = SENS1.m(i);
                  Save_theta_R_D(i,s,nu) = SENS2.m(i);
                  Save_theta_P_W(k,s,nu) = SENS3.m(k);
                  Save_theta_delta_DA(n,s,nu) = SENS4.m(n);

                  Save_C_RT(s,nu) = Z_S.l
         );
         UB(nu) = sum(i, Save_R_U(i,nu)*C_dw(i)) + sum(i, Save_R_D(i,nu)*C_up(i)) + sum(i, Save_P(i,nu)*C(i)) + sum(s, pi(s) * Save_C_RT(s,nu));

         Save_LB(nu,t) = LB(nu);
         Save_UB(nu,t) = UB(nu);
*conv=abs((UB(nu)-LB(nu))/UB(nu))*100;
);
*);
);

Display Save_DD_U, Save_DD_D, Save_LB, Save_UB;

File Upper_bound / Upper_bound.txt /;
put Upper_bound;
loop((nu,t),
  put nu.tl, t.tl, Save_UB(nu,t) /
);
putclose;
File Lower_bound / Lower_bound.txt /;
put Lower_bound;
loop((nu,t),
  put nu.tl, t.tl, Save_LB(nu,t) /
);
putclose;
File Down_reserve / Down_reserve.txt /;
put Down_reserve;
loop((nu,z,t),
  put nu.tl, z.tl, t.tl, Save_DD_D(nu,z,t) /
);
putclose;
File Up_reserve / Up_reserve.txt /;
put Up_reserve;
loop((nu,z,t),
  put nu.tl, z.tl, t.tl, Save_DD_U(nu,z,t) /
);
putclose;
