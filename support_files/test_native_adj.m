function [native_adj_scores, corr_val] = test_native_adj(patient_conn, patient_roi, atlas_conn, atlas_std, all_inds)

    num_elecs = size(patient_conn,1);
    
    native_adj_scores = NaN*zeros(num_elecs);
    
    pred_mat = NaN(num_elecs);
    
    % double loop through electrodes
    for i = 1:num_elecs
        for j = 1:num_elecs
            % get regions from patient roi
            region1 = patient_roi(i);
            region2 = patient_roi(j);
            
            if region1==region2
                      
            else
                %atlas region1 
                atlas_row = find(all_inds==region1);
                atlas_col = find(all_inds==region2);

                atlas_edge_val = atlas_conn(atlas_row, atlas_col);
                atlas_var = atlas_std(atlas_row, atlas_col);

                new_score = (patient_conn(i,j)-atlas_edge_val)./atlas_var;

                if isempty(new_score)
                    native_adj_scores(i,j) = NaN;
                    pred_mat(i,j) = NaN;
                else
                    native_adj_scores(i,j) = new_score;
                    pred_mat(i,j) = atlas_edge_val;
                end
            end
        end
    end
    
    % find predictions which are not NaN
    non_nan_preds = find(~isnan(pred_mat));
    
    % set up real vectors of real connectivity and predicted connectivity
    real_conn = patient_conn(non_nan_preds);
    pred_conn = pred_mat(non_nan_preds);
    
    
    
    % get a global correlation
    corr_val = corr(real_conn, pred_conn);
    
    % get number of edges
    num_real_edges = length(real_conn);
    
end