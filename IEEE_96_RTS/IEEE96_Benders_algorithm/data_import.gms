option decimals=1;
Sets
         i       'set of conventional plants'    /i1*i96/
         j       'set of loads'                  /j1*j3/
         k       'set of wind farms'             /k1*k3/
         n       'set of nodes'                  /n1*n3/
         s       'set of scenarios'              /s1*s100/
         t       'set of time units'             /t1*t24/
         y       'for sos1 linearization'        /y1*y2/
         z       'set of reserve control zones'  /z1*z3/

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
 C_mn/
;
Table datagen(i,dataGENset) 'generators data'
$ondelim
$include data_gen.csv
$offdelim

Table P_b(s,t) 'actual wind generation quantities in scenario s, MWh'
$ondelim
$include data_wind_scenarios.csv
$offdelim

Table D_data(j,t) 'demand of load j, MWh'
$ondelim
$include data_demand.csv
$offdelim

Table ML(j,n) '1/0 matrix for nodal location of load'
$ondelim
$include data_load_location.csv
$offdelim

Parameter X(n,n) 'line reactance, p.u.'
/
$ondelim
$include data_reactances.csv
$offdelim
/;

Parameter T_cap(n,n) 'transmission line capacity, MW'
/
$ondelim
$include data_trans_capacity.csv
$offdelim
/;

Parameter MN(n,n) '1/0 matrix for adjacent nodes'
/
$ondelim
$include data_node_connection.csv
$offdelim
/;

Table MW(k,n) '1/0 matrix for nodal location of wind farms'
$ondelim
$include data_wind_location.csv
$offdelim

Table MG(i,n) '1/0 matrix for nodal location of generators'
$ondelim
$include data_gen_location.csv
$offdelim

Table MZ(i,z) '1/0 matrix for zonal location of generators'
$ondelim
$include data_gen_zones.csv
$offdelim

Parameter wind_cap(k) 'wind capacity'
/
$ondelim
$include data_wind_cap.csv
$offdelim
/;

*$ontext
Table DDUP(t,z) 'reserve requirments data'
$ondelim
$include data_res_UP.csv
$offdelim

Table DDDW(t,z) 'reserve requirments data'
$ondelim
$include data_res_DW.csv
$offdelim
*$offtext

Parameter        P_max(i);
P_max(i)=datagen(i,'P_max');
Parameter        P_min(i);
P_min(i)=datagen(i,'P_min');
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

Parameter P_a(k,t,s)     'actual wind generation quantities in scenario s, MW';
Scalar MM /100000/;
Scalar   VOLL            'value of lost load, Euro/MWh'        /500/;
Scalar   VOSW            'value of spilled wind, Euro/MWh'     /1000/;
Parameter pi(s)          'probability of scenario s';
pi(s) = 1/card(s);
Parameter        W_exp_data(k,t)        'expectation of wind power production';

Display MZ;
