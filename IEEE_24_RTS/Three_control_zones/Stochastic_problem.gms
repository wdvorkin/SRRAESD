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
option iterlim=999999999;     // avoid limit on iterations
*option reslim=300;            // timelimit for solver in sec.
option optcr=0;             // gap tolerance
option solprint=ON;           // include solution print in .lst file
option limrow=100;            // limit number of rows in .lst file
option limcol=100;            // limit number of columns in .lst file
option decimals=1;            // defines decimals in output file
//--------------------------------------------------------------------

$Include data_import.gms

Variables
         z               'expected system costs'
*Day-ahead decision variables
         p(i,t)          'production quantity for power plant i, MWh'
         p_w(k,t)        'production quantity for wind farm k, MWh'
         delta_d(n,t)    'node angle'
         R_U(i,t)        'accepted up-reserve quantity of unit i'
         R_D(i,t)        'accepted down-reserve quantity of unit i'
         DD_RU(t)        'auxiliary up-demand requirements variable'
         DD_RD(t)        'auxiliary up-demand requirements variable'
*Balancing decision variables
         rr_u(i,t,s)       'up-reserve deployment, MWh'
         rr_d(i,t,s)       'down-reserve deployment, MWh'
         w_sp(k,t,s)       'wind spillage, MWh'
         l(j,t,s)          'load shedding, MWh'
         delta_o(n,t,s)    'node angle'
Positive variables       p, p_w, rr_u, rr_d, w_sp, l, R_U, R_D, DD_RU, DD_RD;

*reference bus
delta_d.FX('n1',t)=0;
delta_o.FX('n1',t,s)=0;

Equations
         OBJ_SD          'expected system costs'
*Day-ahead constraints
         URB(t)          'up-reserve balance'
         DRB(t)          'down-reserve balance'
         URLd(i,t)       'up-reserve limit at day-ahead stage'
         DRLd(i,t)       'up-reserve limit at day-ahead stage'
         PB_da(n,t)      'day-ahead power balance'
         MaPC_c(i,t)     'maximum power capacity of convetional power plants'
         MiPC_c(i,t)     'minimum power capacity of convetional power plants'
         MPC_w(k,t)      'maximum power capacity of wind farms'
         MTC_d(n,m,t)    'maximum transmission capacity at day-ahead stage'
*Balancing constraints
         PB_o(n,t,s)             'balancing stage power balance'
         URDL(i,t,s)             'up-reserve deployment limit'
         DRDL(i,t,s)             'down-reserve deployment limit'
         MTC_o(n,m,t,s)          'maximum transmission capacity at balancing stage'
         WSL(k,t,s)              'wind spillage limit for wind farm k in scenario s'
         SHL(j,t,s)              'load shedding limit'
;
OBJ_SD..
         z =e=   sum(t,
                 sum(i, C(i)*p(i,t) + C_up(i)*R_U(i,t) + C_dw(i)*R_D(i,t))) +
                 sum(t,
                 sum(s, pi(s)*(
                 sum(i, C_pl(i)*rr_u(i,t,s)- C_mn(i)*rr_d(i,t,s)) + sum(j, VOLL*l(j,t,s)))))
                 ;
URB(t)..
         sum(i, R_U(i,t)) =e= DD_RU(t);
DRB(t)..
         sum(i, R_D(i,t)) =e= DD_RD(t);
URLd(i,t)..
         R_U(i,t) =l= R_up_max(i);
DRLd(i,t)..
         R_D(i,t) =l= R_dw_max(i);
PB_da(n,t)..
         sum(i$MG(i,n), p(i,t)) + sum(k$MW(k,n), p_w(k,t))  - sum(j$ML(j,n), D(j,t))
         =e=
         sum(m$(MN(n,m)=1),(delta_d(n,t)-delta_d(m,t))/X(n,m));
MaPC_c(i,t)..
         p(i,t) =l= P_max(i) - R_U(i,t);
MiPC_c(i,t)..
         p(i,t) =g= R_D(i,t);
MPC_w(k,t)..
         p_w(k,t) =l= smax(s, P_a(k,t,s));
MTC_d(n,m,t)$MN(n,m)..
         (1/X(n,m))*(delta_d(n,t)-delta_d(m,t)) =l= T_cap(n,m);
PB_o(n,t,s)..
         sum(i$MG(i,n), rr_u(i,t,s)-rr_d(i,t,s)) + sum(j$ML(j,n), l(j,t,s))
       + sum(k$MW(k,n), P_a(k,t,s)) - sum(k$MW(k,n), p_w(k,t))
       - sum(k$MW(k,n), w_sp(k,t,s)) =e=
         sum(m$(MN(n,m)=1),(1/X(n,m))*(delta_o(n,t,s)-delta_d(n,t)-delta_o(m,t,s)+delta_d(m,t)));
URDL(i,t,s)..
         rr_u(i,t,s) =l= R_U(i,t);
DRDL(i,t,s)..
         rr_d(i,t,s) =l= R_D(i,t);
MTC_o(n,m,t,s)$MN(n,m) ..
         (1/X(n,m))*(delta_o(n,t,s)-delta_o(m,t,s)) =l= T_cap(n,m);
WSL(k,t,s)..
         w_sp(k,t,s) =l= P_a(k,t,s);
SHL(j,t,s)..
         l(j,t,s) =l= D(j,t);
Model StochD /all/;

File Stoch_output/Stoch_output.txt/;
put Stoch_output;
Set       alpha   'set of wind penetration levels'    /alpha0/;
Parameter Alpha_increment(alpha);
Parameter Penetration_level(alpha);

Parameter Total_expected_costs_for_diff_alpha(alpha);
Parameter Total_expected_costs_for_diff_alpha_and_t(alpha,t);

Loop(alpha,
         P_a(k,t,s)=P_b(s,t)*100;
         W_exp(k,t)=sum(s, pi(s)*P_a(k,t,s));
         Penetration_level(alpha) = sum(k, smax(s, P_a(k,'t19',s)))/sum(j, D(j,'t19'))*100;

         Solve StochD using lp minimizing z;

         Total_expected_costs_for_diff_alpha(alpha)=z.l + 0.000001;
         Total_expected_costs_for_diff_alpha_and_t(alpha,t) = sum((i), C_up(i)*R_U.l(i,t) + C_dw(i)*R_D.l(i,t)) + sum((i), C(i)*p.l(i,t)) + sum(s, pi(s)*(sum(i, C_pl(i)*rr_u.l(i,t,s)- C_mn(i)*rr_d.l(i,t,s)) + sum(j, VOLL*l.l(j,t,s))))  + 0.000001;

         loop(t,
         put ord(t), Total_expected_costs_for_diff_alpha_and_t(alpha,t)/;
         );
);
Display Total_expected_costs_for_diff_alpha, Total_expected_costs_for_diff_alpha_and_t, DD_RU.l, DD_RD.l;
