clear
close all
sbioreset
clc
%% Step 1. Generate plausible patients and model prediction
% LHS sampling
immune_oncology_model_NSCLC
%%
n_PSA = 15;
params_in  = PSA_param_in_NSCLC;
params_in  = PSA_setup(model,params_in,n_PSA);

trial_name = 'durva';

%% Add PK variability

load('PK_pars.mat')
names.pars = {};
names.pars = [names.pars; 'q_P_aPDL1'; 'k_cl_aPDL1'; 'gamma_P_aPDL1'; 'V_C'; 'k_cln_aPDL1'; 'Kc_aPDL1'];
names.labels = {};
names.labels = [names.labels; 'Capillary filtration rate of aPDL1'; 'Linear clearance rate of aPDL1';...
                              'Volume fraction of interstitium in Vp'; 'Total blood volume';...
                              'Non-linear clearance rate of aPDL1'; 'Michaelis–Menten constant for non-linear CL'];
params_in = PSA_setup_PK(params_in,n_PSA,names,psi,V,opt_avg);

%% Simulate for initial condition (pre-treatment)
SIMDATA = cell(1,n_PSA);
for i = 1:n_PSA
    display(['Sample ',num2str(i),'/',num2str(n_PSA)]);
    % Set the new parameters
    if i > 1
        delete(model_PSA)
    end
    model_PSA = copyobj(model);
    variantObj = addvariant(model_PSA, ['v',num2str(i,'%5.5i')]);
    for j = 1:length(params_in.names)
        if ~isempty(sbioselect (model, 'Type', 'parameter', 'Name', params_in.names{j}))
            addcontent(variantObj, {'parameter', params_in.names{j}, 'Value', params_in.(params_in.names{j}).LHS(i)});
        elseif ~isempty(sbioselect (model, 'Type', 'compartment', 'Name', params_in.names{j}))
            addcontent(variantObj, {'compartment', params_in.names{j}, 'Capacity', params_in.(params_in.names{j}).LHS(i)});
        else
            disp(['Unable to identify parameter/compartment named ', params_in.names{j}, ' in the model'])
        end
    end
    % Set Initial Conditions
    for k = 1:length(variantObj.content)
        variantParam = variantObj.content{k};
        if strcmp(variantParam(2),'initial_tumour_diameter')
            % Get Initial Tumour Diameter from 'variant'
            D_Ti = sbioselect(model,'Name','initial_tumour_diameter');
            tumour_diameter.Value = variantParam{4};
            tumour_diameter.Units = D_Ti.ValueUnits;
            tumour_diameter.Notes = D_Ti.Notes;
        end
    end
    % Calculate Target Tumour Volume
    tumour_volume = 4/3*pi*(tumour_diameter/2)^3;
    % Reset Output Times for IC Simulation (assumes time in days)
    config = getconfigset(model);
    set(config,'StopTime',8000);
    set(config.SolverOptions,'OutputTimes',0:8000);
    set(config.SolverOptions,'MaxStep',1);

    simData = sbiosimulate(model,variantObj);
    SIMDATA{i,1} = simData;
end

%%
vars = {'V_T','V_T.Mac_M1','V_T.Mac_M2','V_T.T1','V_T.T1_exh','V_T.T0','V_T.Th','V_T.Th_exh','V_T.C1'};
for i = 1:n_PSA
    numRows = 3; % Number of rows in the grid
    numCols = 3; % Number of columns in the grid

    % Create figure
    figure;
    simData = SIMDATA{i,1};
    for j = 1:length(vars)
        var = vars{j};
        [t,x,~] = selectbyname(simData, var);
        subplot(numRows, numCols, j); % Create subplot in a grid
        plot(t,x); % Plot the data for the current variable
        title(var); % Add a title to the subplot
    end

end

%% Plots for each compartment
% vars = {'V_T','V_T.Mac_M1','V_T.Mac_M2','V_T.T1','V_T.T1_exh','V_T.T0',...
%     'V_T.Th','V_T.Th_exh','V_T.C1','V_T.IFNg','V_T.TGFb','V_T.APC','V_T.mAPC'};

comps = {model.Compartments.Name};
for k = 1:length(comps)
    comp_name = comps{k};
    comp = sbioselect(model,'Name',comp_name,'Type','Compartment');
    vars = {comp.Species.Name};
    vars = append([comp_name '.'],vars);
    num_vars = length(vars);
    figure;
    numRows = ceil(sqrt(num_vars)); % Number of rows in the grid
    numCols = ceil(sqrt(num_vars)); % Number of columns in the grid
    for i = 1:num_vars
        var = vars{i};
        subplot(numRows, numCols, i); % Create subplot in a grid
        % Create figure
        for j = 1:n_PSA
            simData = SIMDATA{j,1};
            [t,x,~] = selectbyname(simData, var);
            plot(t,x); % Plot the data for the current variable
            hold on
            title(var); % Add a title to the subplot
            % set(gca, 'YScale', 'log')
            % legend(append('Patient ',string(1:n_PSA)))
        end
    end
end
% I should look at how all the synaptic species compare within a single
% patient, across compartments. Same with cells, and unbound cytokines as
% well.

%% dosing
config = getconfigset(model);
set(config,'StopTime',400);
set(config.SolverOptions,'OutputTimes',0:400);

dose_schedule = schedule_dosing({'ipilimumab'});

sbioaccelerate(model, dose_schedule)
tic
[simDataPSA, params_out] = simbio_PSA(model,params_in,params_out,dose_schedule);
toc

%% Plots for each compartment after dosing
sd_cell = {simDataPSA.simData};
% comps = {model.Compartments.Name};
comps = {'V_T'};
for k = 1:length(comps)
    comp_name = comps{k};
    comp = sbioselect(model,'Name',comp_name,'Type','Compartment');
    vars = {comp.Species.Name};
    vars = append([comp_name '.'],vars);
    num_vars = length(vars);
    figure;
    numRows = ceil(sqrt(num_vars)); % Number of rows in the grid
    numCols = ceil(sqrt(num_vars)); % Number of columns in the grid
    for i = 1:num_vars
        var = vars{i};
        subplot(numRows, numCols, i); % Create subplot in a grid
        % Create figure
        for j = 1:n_PSA
            simData = sd_cell{j};
            if ~isempty(simData)
                [t,x,~] = selectbyname(simData, var);
                plot(t,x); % Plot the data for the current variable
            end
            hold on
            title(var); % Add a title to the subplot
        end
    end
end