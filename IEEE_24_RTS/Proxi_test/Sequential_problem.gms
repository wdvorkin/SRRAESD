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
         z               'operating reserve costs'
;
Equations
         OBJ_OR          'total operating reserve costs'
         URB(t)          'up-reserve balance'
         DRB (t)         'down-reserve balance'
         URLd(i,t)       'up-reserve limit at day-ahead stage'
         DRLd(i,t)       'up-reserve limit at day-ahead stage'
         UDB(i,t)        'up- and dow-reserve balance for unit i'
;
OBJ_OR..
         sum((i,t), C_up(i)*R_U(i,t) + C_dw(i)*R_D(i,t)) =e= z;
URB(t)..
         sum(i, R_U(i,t)) =e= D_RU(t);
DRB(t)..
         sum(i, R_D(i,t)) =e= D_RD(t);
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

* UC variables
Binary variables
         u(i,t)            'commitment status of conventional units'
Positive variables
         CSU(i,t)
;
Equations
         OBJ_DA          'day-ahead system costs'
         PB_da(n,t)      'day-ahead power balance'
         MaPC_c(i,t)     'maximum power capacity of convetional power plants'
         MiPC_c(i,t)     'minimum power capacity of convetional power plants'
         MPC_w(k,t)      'maximum power capacity of wind farms'
         MTC(n,m,t)      'maximum transmission capacity'

         MiPC_c2(i,t)    'minimum power capacity of convetional power plants'
         STUP(i,t)       'strat-up binary logic'
         STUP_in(i,t)    'strat-up binary logic'
         RUC_in(i,t)     'ramping up constraint'
         RUC(i,t)        'ramping up constraint'
         RDC_in(i,t)     'ramping down constraint'
         RDC(i,t)        'ramping down constraint'
;
OBJ_DA..
         zz =e= sum((i,t), C(i)*p(i,t) + CSU(i,t));
PB_da(n,t)..
         sum(i$MG(i,n), p(i,t)) + sum(k$MW(k,n), p_w(k,t))  - sum(j$ML(j,n), D(j,t))
         =e=
         sum(m$(MN(n,m)=1),(delta(n,t)-delta(m,t))/X(n,m));
MaPC_c(i,t)..
         p(i,t) =l= [P_max(i) - R_U.l(i,t)]*u(i,t);
MiPC_c(i,t)..
         p(i,t) =g= R_D.l(i,t);
MiPC_c2(i,t)..
         p(i,t) =g= P_min(i)*u(i,t);
MPC_w(k,t)..
         p_w(k,t) =l= W_exp(k,t);
MTC(n,m,t)$MN(n,m)..
         (1/X(n,m))*(delta(n,t)-delta(m,t)) =l= T_cap(n,m);
STUP(i,t)$(ord(t)>1)..
         CSU(i,t) =G= C_su(i) * (u(i,t)-u(i,t-1));
STUP_in(i,t)$(ord(t)=1)..
         CSU(i,t) =G= C_su(i) * (u(i,t)-u_init(i));
RUC(i,t)$(ord(t)>1)..
         p(i,t) - p(i,t-1) =L= Ramp_up(i);
RUC_in(i,t)$(ord(t)=1)..
         p(i,t) - P_init(i) =L= Ramp_up(i);
RDC(i,t)$(ord(t)>1)..
         p(i,t-1) - p(i,t) =L= Ramp_dw(i);
RDC_in(i,t)$(ord(t)=1)..
         P_init(i) - p(i,t) =L= Ramp_dw(i);
Model Conv_DA /OBJ_DA, PB_da, MaPC_c, MiPC_c2, MiPC_c, MPC_w, MTC, STUP, STUP_in, RUC, RUC_in, RDC, RDC_in/;

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
         sum((s,t), pi(s)*(
         sum(i, C_pl(i)*rr_u(i,t,s)- C_mn(i)*rr_d(i,t,s)) + sum(j, VOLL*l(j,t,s))))
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
         l(j,t,s) =l= D(j,t);
Model ConvD_RT /OBJ_RT, PB_o, URDL, DRDL, MTC_o, WSL, SHL/;

P_a(k,t,s) = P_b(s,t)*200;
W_exp(k,t) = sum(s, pi(s)*P_a(k,t,s) );

* Probabilistic requirements
* 90%
*D_RU(t) = min(sum(i, R_up_max(i)) , (sum(k, W_exp(k,t)) - sum(k, P_a(k,t,'s6'))));
*D_RD(t) = min(sum(i, R_dw_max(i))-60, (sum(k, P_a(k,t,'s95')) - sum(k, W_exp(k,t))));

* 95%
D_RU(t) = min(sum(i, R_up_max(i)) , (sum(k, W_exp(k,t)) - sum(k, P_a(k,t,'s3'))));
D_RD(t) = min(sum(i, R_dw_max(i))-60, (sum(k, P_a(k,t,'s98')) - sum(k, W_exp(k,t))));

* Enhanced requirements
*D_RU(t) = datares(t,'D_RU');
*D_RD(t) = datares(t,'D_RD');

Solve Conv_OR using lp minimizing z;
Solve Conv_DA using mip minimizing zz;
Solve ConvD_RT using lp minimizing zzz;

Display u.l, CSU.l, C_su, STUP.m, STUP_in.m, RUC.m, RUC_in.m, RDC.m, RDC_in.m;



*$ontext
Parameter Cost_ts(t);

Cost_ts(t) = sum((i), C_up(i)*R_U.l(i,t) + C_dw(i)*R_D.l(i,t)) + sum((i), C(i)*p.l(i,t) + CSU.l(i,t)) + sum((s), pi(s)*(
sum(i, C_pl(i)*rr_u.l(i,t,s)- C_mn(i)*rr_d.l(i,t,s)) + sum(j, VOLL*l.l(j,t,s))));

Display Cost_ts;


File Results_ts / Results_ts.txt /;
put Results_ts;
loop((t),
  put t.tl, D_RU(t), D_RD(t), Cost_ts(t) /
);
*$offtext





$ontext
*** Wind penetation data
Sets         alpha   'set of wind penetr. levels'    /alpha0/;

Parameter Total_expected_costs_for_diff_alpha(alpha);
Loop(alpha,
         P_a(k,t,s)=P_b(s,t)*100;
         W_exp(k,t)=sum(s, pi(s)*P_a(k,t,s));

         D_RU(t) = min(sum(i, R_up_max(i)) , (sum(k, W_exp(k,t)) - sum(k, P_a(k,t,'s6'))));
         D_RD(t) = min(sum(i, R_dw_max(i))-60 , (sum(k, P_a(k,t,'s95')) - sum(k, W_exp(k,t))));

         Solve Conv_OR using lp minimizing z;
         Solve Conv_DA using lp minimizing zz;
         Solve ConvD_RT using lp minimizing zzz;

         Total_expected_costs_for_diff_alpha(alpha)=z.l + zz.l + zzz.l + 0.000001;
         Display D_RU, D_RD, W_exp,Total_expected_costs_for_diff_alpha;

);

File SeqD_R_U / SeqD_R_U.txt /;
put SeqD_R_U;
loop((i,t),
  put i.tl, t.tl, R_U.l(i,t) /
);
putclose;
File SeqD_R_D / SeqD_R_D.txt /;
put SeqD_R_D;
loop((i,t),
  put i.tl, t.tl, R_D.l(i,t) /
);
putclose;
File SeqD_P / SeqD_P.txt /;
put SeqD_P;
loop((i,t),
  put i.tl, t.tl, p.l(i,t) /
);
putclose;
File SeqD_P_W / SeqD_P_W.txt /;
put SeqD_P_W;
loop((k,t),
  put k.tl, t.tl, p_w.l(k,t) /
);
putclose;
File SeqD_delta/ SeqD_delta.txt /;
put SeqD_delta;
loop((n,t),
  put n.tl, t.tl, delta.l(n,t) /
);
putclose;


Parameter Cong(n,m,t);
Cong(n,m,t) = abs((1/X(n,m))*(delta.l(n,t)-delta.l(m,t)))$(MN(n,m)=1) + 9999999$(MN(n,m) ne 1) ;
Display Cong;
$offtext
