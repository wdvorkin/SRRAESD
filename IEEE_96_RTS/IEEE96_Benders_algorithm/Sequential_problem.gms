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
option iterlim=999999999;           // avoid limit on iterations
option optcr=0;                     // gap tolerance
option solprint=OFF;                 // include solution print in .lst file
option limrow=0;                // limit number of rows in .lst file
option limcol=0;                  // limit number of columns in .lst file
option decimals=1;                  // defines decimals in output file
//--------------------------------------------------------------------

$Include data_import.gms

Positive variables
         R_U(i,t)        'accepted up-reserve quantity of unit i'
         R_D(i,t)        'accepted down-reserve quantity of unit i'
Free variables
         zr              'operating reserve costs'
;
Equations
         OBJ_OR          'total operating reserve costs'
         URB(t,z)          'up-reserve balance'
         DRB(t,z)          'down-reserve balance'
         URLd(i,t)       'up-reserve limit at day-ahead stage'
         DRLd(i,t)       'up-reserve limit at day-ahead stage'
         UDB(i,t)        'up- and dow-reserve balance for unit i'
;

Parameters
D_RU(t,z)
D_RD(t,z)
;

OBJ_OR..
         sum((i,t), C_up(i)*R_U(i,t) + C_dw(i)*R_D(i,t)) =e= zr;
URB(t,z)..
         sum(i$MZ(i,z), R_U(i,t)) =E= D_RU(t,z);
DRB(t,z)..
         sum(i$MZ(i,z), R_D(i,t)) =E= D_RD(t,z);
URLd(i,t)..
         R_U(i,t) =l= R_up_max(i);
DRLd(i,t)..
         R_D(i,t) =l= R_dw_max(i);
UDB(i,t)..
         R_U(i,t) + R_D(i,t) =L= P_max(i);
Model Conv_OR /OBJ_OR, URB, DRB, URLd, DRLd, UDB/;

Variables
         p(i,t)          'production quantity for power plant i, MWh'
         p_w(k,t)        'production quantity for wind farm k, MWh'
         delta(n,t)      'voltage angle'
         zz              'day-ahead system costs'
Positive variables       p, p_w;
delta.FX('n1',t)=0;
Equations
         OBJ_DA          'day-ahead system costs'
         PB_da(n,t)      'day-ahead power balance'
         MaPC_c(i,t)     'maximum power capacity of convetional power plants'
         MiPC_c(i,t)     'minimum power capacity of convetional power plants'
         MPC_w(k,t)      'maximum power capacity of wind farms'
         MTC(n,m,t)      'maximum transmission capacity'
;
OBJ_DA..
         zz =e= sum((i,t), C(i)*p(i,t));
PB_da(n,t)..
         sum(i$MG(i,n), p(i,t)) + sum(k$MW(k,n), p_w(k,t))  - sum(j$ML(j,n), D_data(j,t))
         =e=
         sum(m$(MN(n,m)=1),(delta(n,t)-delta(m,t))/X(n,m));
MaPC_c(i,t)..
         p(i,t) =l= P_max(i) - R_U.l(i,t);
MiPC_c(i,t)..
         p(i,t) =g= R_D.l(i,t);
MPC_w(k,t)..
         p_w(k,t) =l= W_exp_data(k,t);
MTC(n,m,t)$MN(n,m)..
         (1/X(n,m))*(delta(n,t)-delta(m,t)) =l= T_cap(n,m);
Model Conv_DA /OBJ_DA, PB_da, MaPC_c, MiPC_c, MPC_w, MTC/;

Variables
         rr_u(i,t,s)             'up-reserve deployment, MW'
         rr_d(i,t,s)             'down-reserve deployment, MW'
         w_sp(k,t,s)             'wind spillage, MW'
         l(j,t,s)                'load shedding, MW'
         delta_o(n,t,s)          'voltage angle'
         zzz                     'expected real-time system costs'
Positive variables rr_u, rr_d, w_sp, l;
*reference bus
delta_o.FX('n1',t,s)=0;

Equations
         OBJ_RT                  'expected balancing costs'
         PB_o(n,t,s)             'balancing stage power balance'
         URDL(i,t,s)             'up-reserve deployment limit'
         DRDL(i,t,s)             'down-reserve deployment limit'
         MTC_o(n,m,t,s)          'maximum transmission capacity at balancing stage'
         WSL(k,t,s)              'wind spillage limit for wind farm k in scenario s'
         SHL(j,t,s)              'load shedding limit'
;
OBJ_RT..
         sum(t,
         sum(s, pi(s)*(
         sum(i, C_pl(i)*rr_u(i,t,s)- C_mn(i)*rr_d(i,t,s)) + sum(j, VOLL*l(j,t,s)))))
         =e= zzz;
PB_o(n,t,s)..
         sum(i$MG(i,n), rr_u(i,t,s)-rr_d(i,t,s)) + sum(j$ML(j,n), l(j,t,s))
       + sum(k$MW(k,n), P_a(k,t,s)) - sum(k$MW(k,n), p_w.l(k,t))
       - sum(k$MW(k,n), w_sp(k,t,s)) =e=
         sum(m$(MN(n,m)=1),(1/X(n,m))*(delta_o(n,t,s)-delta.l(n,t)-delta_o(m,t,s)+delta.l(m,t)));

URDL(i,t,s)..
         rr_u(i,t,s) =l= R_U.l(i,t);
DRDL(i,t,s)..
         rr_d(i,t,s) =l= R_D.l(i,t);
MTC_o(n,m,t,s)$MN(n,m) ..
         (1/X(n,m))*(delta_o(n,t,s)-delta_o(m,t,s)) =l= T_cap(n,m);
WSL(k,t,s)..
         w_sp(k,t,s) =l= P_a(k,t,s);
SHL(j,t,s)..
         l(j,t,s) =l= D_data(j,t);
Model ConvD_RT /OBJ_RT, PB_o, URDL, DRDL, MTC_o, WSL, SHL/;

*** Wind penetation data
Sets         alpha   'set of wind penetr. levels'    /alpha0/;

Parameter Total_expected_costs_for_diff_alpha(alpha);
Loop(alpha,
         P_a(k,t,s)=P_b(s,t)*wind_cap(k)*0.25;
         W_exp_data(k,t)=sum(s, pi(s)*P_a(k,t,s));

$ontext
D_RU(t,'z1') = min(W_exp_data('k1',t) - P_a('k1',t,'s1'),  sum(i$MZ(i,'z1'), R_up_max(i)));
D_RU(t,'z2') = min(W_exp_data('k2',t) - P_a('k2',t,'s1'),  sum(i$MZ(i,'z2'), R_up_max(i)));
D_RU(t,'z3') = min(W_exp_data('k3',t) - P_a('k3',t,'s1'),  sum(i$MZ(i,'z3'), R_up_max(i)));
D_RD(t,'z1') = min(- W_exp_data('k1',t) + P_a('k1',t,'s99'),   sum(i$MZ(i,'z1'), R_dw_max(i)));
D_RD(t,'z2') = min(- W_exp_data('k2',t) + P_a('k2',t,'s99'),   sum(i$MZ(i,'z2'), R_dw_max(i)));
D_RD(t,'z3') = min(- W_exp_data('k3',t) + P_a('k3',t,'s99'),   sum(i$MZ(i,'z3'), R_dw_max(i)));
$offtext

$ontext
D_RU(t,'z1') = min(W_exp_data('k1',t) - P_a('k1',t,'s6'),  sum(i$MZ(i,'z1'), R_up_max(i)));
D_RU(t,'z2') = min(W_exp_data('k2',t) - P_a('k2',t,'s6'),  sum(i$MZ(i,'z2'), R_up_max(i)));
D_RU(t,'z3') = min(W_exp_data('k3',t) - P_a('k3',t,'s6'),  sum(i$MZ(i,'z3'), R_up_max(i)));
D_RD(t,'z1') = min(- W_exp_data('k1',t) + P_a('k1',t,'s95'),   sum(i$MZ(i,'z1'), R_dw_max(i)));
D_RD(t,'z2') = min(- W_exp_data('k2',t) + P_a('k2',t,'s95'),   sum(i$MZ(i,'z2'), R_dw_max(i)));
D_RD(t,'z3') = min(- W_exp_data('k3',t) + P_a('k3',t,'s95'),   sum(i$MZ(i,'z3'), R_dw_max(i)));
$offtext

*$ontext
D_RU(t,z) = DDUP(t,z);
D_RD(t,z) = DDDW(t,z);
*$offtext

         Solve Conv_OR using lp minimizing zr;
         Solve Conv_DA using lp minimizing zz;
         Solve ConvD_RT using lp minimizing zzz;

         Total_expected_costs_for_diff_alpha(alpha)=zr.l + zz.l + zzz.l;
         Display zr.l, zz.l, zzz.l, R_U.l, R_D.l;
);

Parameter total_costs;
total_costs = zr.l + zz.l + zzz.l;

Display total_costs;

Parameter total_cost_t(t);
total_cost_t(t) = sum((i), C_up(i)*R_U.l(i,t) + C_dw(i)*R_D.l(i,t)) + sum((i), C(i)*p.l(i,t)) + sum(s, pi(s)*( sum(i, C_pl(i)*rr_u.l(i,t,s)- C_mn(i)*rr_d.l(i,t,s)) + sum(j, VOLL*l.l(j,t,s))));

Display total_cost_t;
Display D_RU, D_RD;
