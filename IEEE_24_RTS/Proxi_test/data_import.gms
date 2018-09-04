option decimals=1;
Sets
         i       'set of conventional plants'    /i1*i12/
         j       'set of loads'                  /j1*j17/
         k       'set of wind farms'             /k1*k6/
         n       'set of nodes'                  /n1*n24/
         s       'set of scenarios'              /s1*s100/
         t       'set of time units'             /t1*t24/
         y       'for sos1 linearization'        /y1*y2/

Alias(m,n)
Sets
dataGENset  'set of import generation data parameters'
/P_max
P_min
C
C_up
C_dw
R_up_max
R_dw_max
C_pl
C_mn
Ramp_up
Ramp_dw
P_init
C_su
u_init/
dataRESset
/D_RU
 D_RD/
;
Parameters
D_RU(t)         'up-reserve demand'
D_RD(t)         'down-reserve demand'

Table datagen(i,dataGENset) 'generators data'
$ondelim
$include data_gen.csv
$offdelim

Table P_b(s,t) 'actual wind generation quantities in scenario s, MWh'
$ondelim
$include data_wind_scenarios.csv
$offdelim

Table D(j,t) 'demand of load j, MWh'
$ondelim
$include data_demand.csv
$offdelim

Table datares(t,dataRESset) 'reserve requirments data'
$ondelim
$include data_res.csv
$offdelim

Table ML(j,n) '1/0 matrix for nodal location of load'
$ondelim
$include data_load_location.csv
$offdelim

Table X(n,n) 'line reactance, p.u.'
$ondelim
$include data_reactances.csv
$offdelim

Table T_cap(n,n) 'transmission line capacity, MW'
$ondelim
$include data_trans_capacity.csv
$offdelim

Table MN(n,n) '1/0 matrix for adjacent nodes'
$ondelim
$include data_node_connection.csv
$offdelim

Table MW(k,n) '1/0 matrix for nodal location of wind farms'
$ondelim
$include data_wind_location.csv
$offdelim

Table MG(i,n) '1/0 matrix for nodal location of generators'
$ondelim
$include data_gen_location.csv
$offdelim

Parameters
P_max(i)
P_min(i)
C(i)
C_up(i)
C_dw(i)
R_up_max(i)
R_dw_max(i)
C_pl(i)
C_mn(i)
Ramp_up(i)
Ramp_dw(i)
P_init(i)
C_su(i)
u_init(i)
;
P_max(i)=datagen(i,'P_max');
P_min(i)=datagen(i,'P_min');
C(i)=datagen(i,'C');
C_up(i)=datagen(i,'C_up');
C_dw(i)=datagen(i,'C_dw');
R_up_max(i)=datagen(i,'R_up_max');
R_dw_max(i)=datagen(i,'R_dw_max');
C_pl(i)=datagen(i,'C_pl');
C_mn(i)=datagen(i,'C_mn');
Ramp_up(i)=datagen(i,'Ramp_up');
Ramp_dw(i)=datagen(i,'Ramp_dw');
P_init(i)=datagen(i,'P_init');
C_su(i)=datagen(i,'C_su');
u_init(i)=datagen(i,'u_init');

Parameter P_a(k,t,s)     'actual wind generation quantities in scenario s, MW';
Scalar   VOLL             'value of lost load, Euro/MWh'  /500/;
Scalar   VOSW             'value of spilled wind, Euro/MWh'  /1000/;
Parameter pi(s)          'probability of scenario s';
pi(s) = 1/card(s);
Parameter        W_exp(k,t)        'expectation of wind power production';
