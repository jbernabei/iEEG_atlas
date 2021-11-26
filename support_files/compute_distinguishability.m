function [D_rs_vec] = compute_distinguishability(adj_matrices,resected_elecs,metric_type,WM_eliminate, patient_roi)

    % calculate number of patients
    num_patients = length(adj_matrices);
    
    % loop through patients
    for pt = 1:num_patients
        
        patient_WM_inds = find(patient_roi{pt}==9171);
        
        % extract resected electrode indices
        pt_res_elecs = resected_elecs{pt};
        
        % loop through frequency bands
        for f = 1:5
            
            % extract adjacency matrix
            patient_adj = adj_matrices{pt}(f).data;
            resected_elec_bool = zeros(size(patient_adj,1),1);
            if WM_eliminate
                patient_adj(patient_WM_inds,:) = [];
                patient_adj(:,patient_WM_inds) = [];
                            % do resected elecs
                
                resected_elec_bool(resected_elecs{pt}) = 1;
                resected_elec_bool(patient_WM_inds) = [];
                pt_res_elecs = find(resected_elec_bool);
            end

            if strcmp(metric_type,'node_strength')
                pt_metric{pt}.data(:,f) = zscore(sum(patient_adj));
            elseif strcmp(metric_type,'betweenness_centrality')
                pt_metric{pt}.data(:,f) = zscore(betweenness_wei(patient_adj));
            elseif strcmp(metric_type,'control_centrality')
                pt_metric{pt}.data(:,f) = zscore(control_centrality(patient_adj));
            elseif strcmp(metric_type,'participation_coefficient')
                [S,Q] = modularity_und(patient_adj,1);
                pt_metric{pt}.data(:,f) = participation_coef(patient_adj,S);
            end
            
            res_result = pt_metric{pt}.data(pt_res_elecs,f);
            non_res_result = pt_metric{pt}.data(:,f);
            non_res_result(pt_res_elecs) = [];
            
            mwu_result = mwwtest(res_result',non_res_result');
            
            D_rs_vec(pt,f) = 1-max(mwu_result.U(2))./(length(res_result).*length(non_res_result));
            
        end
    end

end