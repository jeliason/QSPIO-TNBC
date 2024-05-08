function [success,idx,M1,M2,Teff,T1exh,T0,Th,Thexh,C_total,VDT] = PSA_NSCLC_iteration(param_struct)
% function [M1_M2,Treg_CD8,CD8_CD4] = PSA_NSCLC_iteration(param_struct)

    load('PSA_iteration.mat','model');
    dt = 56;
    % ub = [6 2.3 2.0 50]'; lb = [.8 0 0 .06]';
    % Set the new parameters
    disp('starting sampling');
    if i > 1
        delete(model_PSA)
    end
    model_PSA = copyobj(model);
    variantObj = addvariant(model_PSA, 'variant');
    fields = fieldnames(param_struct);

    for j = 1:numel(fields)
        param_name = fields{j};
        param_value = param_struct.(param_name);
        if ~isempty(sbioselect (model, 'Type', 'parameter', 'Name', param_name))
            addcontent(variantObj, {'parameter', param_name, 'Value', param_value});
        elseif ~isempty(sbioselect (model, 'Type', 'compartment', 'Name', param_name))
            addcontent(variantObj, {'compartment', param_name, 'Capacity', param_value});
        else
            disp(['Unable to identify parameter/compartment named ', params_name, ' in the model'])
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
    
    try
        % IC Simulation
        simData = sbiosimulate(model,variantObj);
        % Get Tumour Volume from Simulation
        [~,V_T,~] = selectbyname(simData, 'V_T'); % milliliter
        Vt = V_T';
        target_V_T = tumour_volume.Value; % tumour volume in mL
        % Get Time of Target Tumour Size
        idx = find((target_V_T-V_T)<0,1);
        if isempty(idx) % tumour did not reach target size
            disp(['Initial conditions: ' num2str(target_V_T) 'mL not reached (max: ' num2str(max(V_T)) 'mL)']);
            % M1_M2 = inf;
            % Treg_CD8 = inf;
            % CD8_CD4 = inf;
        else
            success = 1;
            [~,M1,~] = selectbyname(simData, 'V_T.Mac_M1');
            [~,M2,~] = selectbyname(simData, 'V_T.Mac_M2');
            [~,Teff,~] = selectbyname(simData, 'V_T.T1');
            [~,T1exh,~] = selectbyname(simData, 'V_T.T1_exh');
            [~,T0,~] = selectbyname(simData, 'V_T.T0');
            [~,Th,~] = selectbyname(simData, 'V_T.Th');
            [~,Thexh,~] = selectbyname(simData, 'V_T.Th_exh');
            [~,C_total,~] = selectbyname(simData, 'V_T.C1');
            CD8 = Teff + T1exh;
            CD4 = T0 + Th + Thexh;
            Treg = T0;
            % Vt = V_T';
            % C_diag = C_total(idx);
            % C_time = C_total;
            M1_M2 = M1(idx)./M2(idx);
            Treg_CD8 = Treg(idx)./CD8(idx);
            CD8_CD4 = CD8(idx)./CD4(idx);
            
            t1 = find(V_T>=target_V_T,1); t2 = t1+dt;
            V1 = V_T(t1); V2 = V_T(t2);
            VDT = dt.*log(2)./(log(V2)-log(V1));
            % filter by lower and upper boundaries of selected model variables
            %             if output(i,1)>ub(1) || output(i,1)<lb(1) || output(i,2)>ub(2) || output(i,2)<lb(2) ||...
            %                     output(i,3)>ub(3) || output(i,3)<lb(3)
            %                 success(i) = 0;
            %                 disp('Does not fall within physiologically reasonable ranges');
            %                 disp(output(i,:))
            %             end
        end
    catch
        disp('There was an error while finding the initial conditions');
    end
end