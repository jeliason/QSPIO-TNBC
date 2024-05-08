function [param_names,params] = create_params(n_PSA)

	%% Step 1. Generate plausible patients and model prediction
	% LHS sampling
	% immune_oncology_model_NSCLC
	load('PSA_iteration1.mat','model');
	% n_PSA = 20;
	% n_PSA = int(n_PSA);

	params_in  = PSA_param_in_NSCLC;
	params_in  = PSA_setup(model,params_in,n_PSA);


	%% Add PK variability

	load('PK_pars.mat','psi','V','opt_avg')
	names.pars = {};
	names.pars = [names.pars; 'q_P_aPDL1'; 'k_cl_aPDL1'; 'gamma_P_aPDL1'; 'V_C'; 'k_cln_aPDL1'; 'Kc_aPDL1'];
	names.labels = {};
	names.labels = [names.labels; 'Capillary filtration rate of aPDL1'; 'Linear clearance rate of aPDL1';...
								'Volume fraction of interstitium in Vp'; 'Total blood volume';...
								'Non-linear clearance rate of aPDL1'; 'Michaelis–Menten constant for non-linear CL'];
	params_in = PSA_setup_PK(params_in,n_PSA,names,psi,V,opt_avg);

	param_names = params_in.names;
	params = params_in.all;
end


