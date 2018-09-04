$ontext
Setting Reserve Requirements to Approximate The Stochastic Dispatch Efficiency
Copyright (C) 2017 Vladimir Dvorkin

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
*option iterlim=999999999;      // avoid limit on iterations
option reslim=100000;           // timelimit for solver in sec.
option optcr=0;                 // gap tolerance
option solprint=OFF;            // include solution print in .lst file
option limrow=100;              // limit number of rows in .lst file
option limcol=100;              // limit number of columns in .lst file
option decimals=1;              // defines decimals in output file
option mip=cplex;

Sets
         i       'set of conventional plants'    /i1*i12/
         j       'set of loads'                  /j1*j17/
         k       'set of wind farms'             /k1*k6/
         n       'set of nodes'                  /n1*n24/
         s       'set of scenarios'              /s1*s100/
         t       'set of time units'             /t1*t24/
         y       'for sos1 linearization'        /y1*y2/
         o       'set of reserve control areas'  /o1*o3/
Alias(m,n)

Parameters
P_b(s)
D(j)
P_a(k,s)
W_exp(k)

Sets
dataGENset  'set of import generation data parameters'
/P_max
 C
 C_up
 C_dw
 R_up_max
 R_dw_max
 C_pl
 C_mn/

dataRESset
/D_RU
 D_RD/
;

Parameters
P_max(i)        'maximum generation capacity, MWh'
C(i)            'marginal production costs of conventional plants, Euro/MWh'
C_up(i)         'up-reserve costs, Euro/MW'
C_dw(i)         'down-reserve costs, Euro/MW'
C_pl(i)         'up-reserve deployment costs, Euro/MW'
C_mn(i)         'down-reserve deployment costs, Euro/MW'
R_up_max(i)     'maximum up-reserve capacity, MW'
R_dw_max(i)     'maximum down-reserve capacity, MW'
D_RU(t)         'up-reserve demand'
D_RD(t)         'down-reserve demand'

Table datagen(i,dataGENset) 'generators data'
$ondelim
$include data_gen.csv
$offdelim
*Display datagen

Table P_b_data(s,t) 'actual wind generation quantities in scenario s, MWh'
$ondelim
$include data_wind_scenarios.csv
$offdelim
*Display P_b_data

Table D_data(j,t) 'demand of load j, MWh'
$ondelim
$include data_demand.csv
$offdelim
*Display D_data

Table ML(j,n) '1/0 matrix for nodal location of load'
$ondelim
$include data_load_location.csv
$offdelim
*Display ML

Table X(n,n) 'line reactance, p.u.'
$ondelim
$include data_reactances.csv
$offdelim
*Display X

Table T_cap(n,n) 'transmission line capacity, MW'
$ondelim
$include data_trans_capacity.csv
$offdelim
*Display T_cap

Table MN(n,n) '1/0 matrix for adjacent nodes'
$ondelim
$include data_node_connection.csv
$offdelim
*Display MN

Table MW(k,n) '1/0 matrix for nodal location of wind farms'
$ondelim
$include data_wind_location.csv
$offdelim
*Display MW

Table MG(i,n) '1/0 matrix for nodal location of generators'
$ondelim
$include data_gen_location.csv
$offdelim
*Display MG

Table MR(i,o) '1/0 matrix for zonal location of generators'
$ondelim
$include data_gen_res_map.csv
$offdelim
Display MR;

Parameter        P_max(i);
P_max(i)=datagen(i,'P_max');
Parameter        C(i);
C(i)=datagen(i,'C');
Parameter        C_up(i);
C_up(i)=datagen(i,'C_up');
Parameter        C_dw(i);
C_dw(i)=datagen(i,'C_dw');
Parameter        R_up_max(i);
R_up_max(i)=datagen(i,'R_up_max');
Parameter        R_dw_max(i);
R_dw_max(i)=datagen(i,'R_dw_max');
Parameter        C_pl(i);
C_pl(i)=datagen(i,'C_pl');
Parameter        C_mn(i);
C_mn(i)=datagen(i,'C_mn');

Parameter P_a_data(k,t,s)     'actual wind generation quantities in scenario s, MW';

Scalar   VOLL             'value of lost load, Euro/MWh'         /500/;
Scalar   VOSW             'value of spilled wind, Euro/MWh'      /1000/;
Parameter pi(s)          'probability of scenario s';
pi(s) = 1/card(s);
Parameter        W_exp_data(k,t)        'expectation of wind power production';
Scalar   MM              'Very huge value'               /5000/;

Variables
         z               'expected system costs'
*Day-ahead decision variables
         p(i)          'production quantity for power plant i, MWh'
         p_w(k)        'production quantity for wind farm k, MWh'
         delta_d(n)    'node angle'
         R_U(i)        'accepted up-reserve quantity of unit i'
         R_D(i)        'accepted down-reserve quantity of unit i'
*Balancing decision variables
         rr_u(i,s)     'up-reserve deployment, MWh'
         rr_d(i,s)     'down-reserve deployment, MWh'
         w_sp(k,s)     'wind spillage, MWh'
         l(j,s)        'load shedding, MWh'
         delta_o(n,s)  'node angle'

Positive variables       p, p_w, R_U, R_D, rr_u, rr_d, w_sp, l;

delta_d.FX('n1')=0;
delta_o.FX('n1',s)=0;

Positive variables
         DD_RU(o)        'maximum up-reserve demand'
         DD_RD(o)        'maximum down-reserve demand'
Free variables
         lambda_RU(o)    'up-reserve price'
         lambda_RD(o)    'down-reserve price'
         lambda_D(n)     'day-ahead price'
Positive variables
         rho_U(i)
         rho_D(i)
         eta_U(i)
         eta_D(i)
         gamma(i)
         mu_U(i)
         mu_L(i)
         phi_U(k)
         phi_L(k)
         mu_delta(n,m)
Binary variables
         u_U(i)
         u_D(i)
         u_U_plus(i)
         u_D_plus(i)
         u_UD(i)
         u_under(i)
         u_over(i)
         u_w_under(k)
         u_w_over(k)
         u_delta(n,m)
;
Equations
         OBJ_BL          'expected system costs'
         PB_o(n,s)     'balancing stage power balance'
         URDL(i,s)     'up-reserve deployment limit'
         DRDL(i,s)     'down-reserve deployment limit'
         MTC_o(n,m,s)  'maximum transmission capacity at balancing stage'
         WSL(k,s)      'wind spillage limit for wind farm k in scenario s'
         SHL(j,s)      'load shedding limit'
         STAT1(i)
         STAT2(i)
         STAT3(i)
         STAT4(k)
         STAT5(n)
         URB(o)          'up-reserve balance'
         DRB(o)          'down-reserve balance'
         URLd(i)       'up-reserve limit at day-ahead stage'
         DRLd(i)       'up-reserve limit at day-ahead stage'
         UDB(i)        'up- and dow-reserve balance for unit i'
         PB_da(n)      'day-ahead power balance'
         MaPC_c(i)     'maximum power capacity of convetional power plants'
         MiPC_c(i)     'minimum power capacity of convetional power plants'
         MPC_w(k)      'maximum power capacity of wind farms'
         MTC(n,m)      'maximum transmission capacity'
         CS11(i)
         CS12(i)
         CS21(i)
         CS22(i)
         CS31(i)
         CS32(i)
         CS41(i)
         CS42(i)
         CS101(i)
         CS102(i)
         CS51(i)
         CS52(i)
         CS61(i)
         CS62(i)
         CS71(k)
         CS72(k)
         CS81(k)
         CS82(k)
         CS91(n,m)
         CS92(n,m)
;
OBJ_BL..
         z =e=   sum(i, C(i)*p(i) + C_up(i)*R_U(i) + C_dw(i)*R_D(i)) +
                 sum(s, pi(s)*(
                 sum(i, C_pl(i)*rr_u(i,s)- C_mn(i)*rr_d(i,s)) + sum(j, VOLL*l(j,s))))
;
PB_o(n,s)..
         sum(i$MG(i,n), rr_u(i,s)-rr_d(i,s)) + sum(j$ML(j,n), l(j,s))
       + sum(k$MW(k,n), P_a(k,s)) - sum(k$MW(k,n), p_w(k))
       - sum(k$MW(k,n), w_sp(k,s)) =e=
         sum(m$(MN(n,m)=1),(1/X(n,m))*(delta_o(n,s)-delta_d(n)-delta_o(m,s)+delta_d(m)));
URDL(i,s)..
         rr_u(i,s) =l= R_U(i);
DRDL(i,s)..
         rr_d(i,s) =l= R_D(i);
MTC_o(n,m,s)$MN(n,m) ..
         (1/X(n,m))*(delta_o(n,s)-delta_o(m,s)) =l= T_cap(n,m);
WSL(k,s)..
         w_sp(k,s) =l= P_a(k,s);
SHL(j,s)..
         l(j,s) =l= D(j);
STAT1(i)..
         C_up(i) - sum(o$MR(i,o), lambda_RU(o)) + gamma(i) + rho_U(i) - eta_U(i) =e= 0;
STAT2(i)..
         C_dw(i) - sum(o$MR(i,o), lambda_RD(o)) + gamma(i) + rho_D(i) - eta_D(i) =e= 0;
STAT3(i)..
         c(i) + mu_U(i) - mu_L(i) - sum(n$MG(i,n), lambda_D(n)) =E= 0;
STAT4(k)..
         phi_U(k) - phi_L(k) - sum(n$MW(k,n), lambda_D(n)) =E= 0;
STAT5(n)..
         sum(m$(MN(n,m)=1), (mu_delta(n,m)+lambda_D(n))/X(n,m) - (mu_delta(m,n)+lambda_D(m))/X(m,n)) =E= 0;
URB(o)..
         sum(i$MR(i,o), R_U(i)) =e= DD_RU(o);
DRB(o)..
         sum(i$MR(i,o), R_D(i)) =e= DD_RD(o);
URLd(i)..
         R_U(i) =l= R_up_max(i);
DRLd(i)..
         R_D(i) =l= R_dw_max(i);
UDB(i)..
         R_U(i) + R_D(i) =L= P_max(i);
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
CS11(i)..
         rho_U(i) =l= MM * u_U(i);
CS12(i)..
         R_up_max(i) - R_U(i) =l= MM * (1-u_U(i));
CS21(i)..
         rho_D(i) =l= MM * u_D(i);
CS22(i)..
         R_dw_max(i) - R_D(i) =l= MM * (1-u_D(i));
CS31(i)..
         eta_U(i) =l= MM * u_U_plus(i);
CS32(i)..
         R_U(i) =l= MM * (1-u_U_plus(i));
CS41(i)..
         eta_D(i) =l= MM * u_D_plus(i);
CS42(i)..
         R_D(i) =l= MM * (1-u_D_plus(i));
CS101(i)..
         gamma(i) =L= MM * u_UD(i);
CS102(i)..
         P_max(i) - R_U(i)- R_D(i) =L= MM *(1- u_UD(i));
CS51(i)..
         mu_U(i) =L= MM * u_over(i);
CS52(i)..
         P_max(i)- R_U(i) - p(i) =L= MM * (1-u_over(i));
CS61(i)..
         mu_L(i) =L= MM * u_under(i);
CS62(i)..
         - R_D(i) + p(i) =L= MM*(1-u_under(i));
CS71(k)..
         phi_U(k) =L= MM * u_w_over(k);
CS72(k)..
         W_exp(k) - p_w(k) =L= MM * (1-u_w_over(k));
CS81(k)..
         phi_L(k) =L= MM * u_w_under(k);
CS82(k)..
         p_w(k) =L= MM * (1-u_w_under(k));
CS91(n,m)$(MN(n,m)=1)..
         mu_delta(n,m) =L= MM * u_delta(n,m);
CS92(n,m)$(MN(n,m)=1)..
         T_cap(n,m) - (1/X(n,m))*(delta_d(n)-delta_d(m)) =L= MM * (1 - u_delta(n,m));
Model    Bilevel /all/;

Option MIP = CPLEX;
$onecho > cplex.opt
epint 1e-15
threads 3
names no
memoryemphasis 1
nodesel 2
$offecho
Bilevel.OptFile = 1;

File Bilevel_output/Bilevel_output.txt/;
put Bilevel_output;

Set       alpha   'set of wind penetr. levels'    /alpha0/;
Parameter Total_expected_costs_for_diff_alpha_and_t(alpha,t);

Loop(t,
         Loop(alpha,
                 P_a_data(k,t,s)=P_b_data(s,t);
                 D(j)=D_data(j,t);
                 P_a(k,s)=P_a_data(k,t,s)*100;
                 W_exp(k)=sum(s, pi(s)*P_a(k,s));

                 Solve    Bilevel using mip minimizing z;
                 put ord(t), DD_RU.l('o1'), DD_RU.l('o2'), DD_RU.l('o3'), DD_RD.l('o1'), DD_RD.l('o2'), DD_RD.l('o3'), z.l/;

                 Total_expected_costs_for_diff_alpha_and_t(alpha,t)=z.l;
         );
);
Display W_exp,Total_expected_costs_for_diff_alpha_and_t, DD_RU.l, DD_RD.l;
