import matlab.engine
import pickle as pkl


eng = matlab.engine.start_matlab()
eng.addpath('scripts/')

param_struct = {
    'k_C1_growth': 0.0160,
    'k_C1_death': 2.0091e-04,
    'k_T1': 0.1108,
    'k_Treg': 0.2179,
    'k_cell_clear': 0.0139,
    'k_C_T1': 1.0133,
    'k_P1_d1': 2.1137e-08,
    'n_T1_clones': 67.0381,
    'n_T0_clones': 296.1959,
    'N_IL2_CD8': 9.7700,
    'N_IL2_CD4': 6.4747,
    'q_T0_T_in': 1.1852e-04,
    'q_T1_T_in': 2.0942e-04,
    'initial_tumour_diameter': 4.2110,
    'T_PD1_total': 1.6162e+04,
    'PD1_50': 24.8685,
    'r_PDL2C1': 0.0011,
    'r_PDL2APC': 0.0307,
    'k_Th_Treg': 0.0043,
    'k_MDSC_rec': 5.0141e+03,
    'k_Mac_rec': 2.6855e+05,
    'k_M1_pol': 0.0149,
   'k_M2_pol': 0.4075,
    'k_CCL2_sec': 2.1667e-12,
    'C_CD47': 63.8422,
    'M_SIRPa': 120.2605,
    'M_PD1_total': 6.5621e+03,
    'k_M1_phago': 0.2636,
    'SIRPa_50': 30.4713,
    'K_Mac_C': 1.0584,
    'q_P_aPDL1': 5.8704e-06,
    'k_cl_aPDL1': 0.3597,
    'gamma_P_aPDL1': 0.0576,
    'V_C': 6.3577,
    'k_cln_aPDL1': 5.1502,
    'Kc_aPDL1': 1.2682
}

success,idx,M1,M2,Teff,T1exh,T0,Th,Thexh,C_total,VDT = eng.PSA_NSCLC_iteration(param_struct,nargout=11)
pkl.dump([success,idx,M1,M2,Teff,T1exh,T0,Th,Thexh,C_total,VDT],open('PSA_NSCLC_iteration.pkl','wb'))

# M1_M2,Treg_CD8,CD8_CD4 = eng.PSA_NSCLC_iteration(param_struct,nargout=3)
