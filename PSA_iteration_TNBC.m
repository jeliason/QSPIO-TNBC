function PSA_iteration_TNBC(start_index, stop_index,save_folder)

immune_oncology_model_TNBC

load('VP.mat','params_in')

%%
config = getconfigset(model);
time = get(config.SolverOptions,'OutputTimes');
n_PSA = length(params_in.(params_in.names{1}).LHS);
model_PSA = [];

for i = start_index:stop_index
    display(['Sample ',num2str(i),'/',num2str(n_PSA)]);
    simDataPSA(i).index = i;
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

    [model_PSA,success,~] = initial_conditions(model_PSA,'Variant',variantObj);
    simDataPSA(i).success = double(success);

    % Run the model with drugs
    if (success)
          try
              simData = sbiosimulate(model_PSA,[],variantObj,[]);
              % Remove simulations that reached the time limit
              if size(simData.Data,1) < length(time)
                  simData = [];
                  simDataPSA(i).success = -1;
                  disp('Simulation takes longer than the preset time limit');
              end
          catch
              disp('Integration tolerance not met');
              simData = [];
              simDataPSA(i).success = -1;
          end

    else
        disp('Initial conditions not reached');
        simData = [];
    end

    % save model output struture within an array of structures
    simDataPSA(i).simData = simData;

    % save model ICs
    simDataPSA(i).ICs = get_ICs(simData);
end
save(save_folder + "/" + string(start_index) + "_" + string(stop_index) + ".mat","simDataPSA");
end
