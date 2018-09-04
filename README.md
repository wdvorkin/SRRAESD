# Online Appendix to the paper "Setting Reserve Requirements to Approximate the Efficiency of the Stochastic Dispatch"

This repository contains the supplementary material for the paper "Setting Reserve Requirements to Approximate The Stochastic Dispatch Efficiency." This is a joint work by Vladimir Dvorkin with the Technical University of Denmark, Stefanos Delikaraoglou with the ETH Zurich, and Juan M. Morales with the University of Malaga.

This repository provides implementations of the sequential dispatch procedure, stochastic dispatch procedure, and proposed bilevel model in the General Algebraic Modeling System (GAMS) software, as well as the data for all simulation instances. This data is provided for the replicability of the proposed models only.  

This repository includes two directories. The directory "IEEE_24_RTS" consists of the data and GAMS code corresponding to experiments in the case study sections V.A--V.E. The directory "IEEE_96_RTS" includes the data and GAMS codes related to the case study section  V.F. The GAMS codes are provided in .gms extension and include the implementation of model (1)-(3) [Sequential_problem.gms], implementation of model (4) [Stochastic_problem.gms], implementation of model (5) [Bilevel_problem.gms], implementation of model (6)-(7) [Benders_algorithm.gms], and data managing script [data_import.gms].
