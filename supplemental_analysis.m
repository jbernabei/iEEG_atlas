% supplmental figures


%% do renderings of original and new atlas
random_colors = rand([40,3]);

color_matrix = ones(117,3);

for i = 1:40
    this_region = custom_atlas{i,2:4};
    for r = 1:length(rmmissing(this_region))
        color_matrix(this_region(r),:) = random_colors(i,:);
    end
end

% writes a file that can be used in BrainNet viewer to visualize 
dlmwrite('atlas_render_color.txt',color_matrix(1:90,:),'delimiter',' ','precision',5)

figure(1);clf;
imagesc(reshape(random_colors,40,1,3))
set(gca,'TickLabelInterpreter', 'none');
yticks([1:40])
yticklabels(custom_atlas{:,1})

%% correlate connectivity between MNI and HUP
for feat = 1:5
        a = conn_band{feat};
        aa = 'coherence';
    
    conn_edge_mni = this_feat_conn(feat).mni;
    conn_edge_hup = this_feat_conn(feat).hup;
    
conn_edge_all = conn_edge_mni(:)+conn_edge_hup(:);
conn_edge_mni(isnan(conn_edge_all)) = [];
conn_edge_hup(isnan(conn_edge_all)) = [];

std_edge_all = std_edge_mni(:)+std_edge_hup(:);
std_edge_mni(isnan(std_edge_all)) = [];
std_edge_hup(isnan(std_edge_all)) = [];

figure(1);
subplot(1,5,feat)
hold on
plot((conn_edge_mni),(conn_edge_hup),'k.','Color',[128 128 128]/255)
[r,p1] = corr((conn_edge_mni)',(conn_edge_hup)')
rvals(feat) = r;
xlabel('MNI connectivity')
ylabel('HUP connectivity')
p = polyfit((conn_edge_mni),(conn_edge_hup),1);
f = polyval(p,(conn_edge_mni));
plot((conn_edge_mni),f,'k-','LineWidth',2)
title(sprintf('%s %s, r = %f, p = %f',a,aa,r,p1))
hold off
end

% figure(2);clf
% bar(rvals)
%% all-node abnormality score all regions
for ch = 1:4717
    if all_possible_labels(ch,2)==0
        grouping{ch,1} = 'uninvolved';
        class_ind(ch,1) = 1;
    elseif all_possible_labels(ch,3)==1
        grouping{ch,1} = 'EIZ';
        class_ind(ch,1) = 2;
    elseif all_possible_labels(ch,4)==1
        grouping{ch,1} = 'SOZ';
        class_ind(ch,1) = 3;
    end
end

new_roi_2 = new_roi;
new_roi_2(incomplete_channels) = [];

single_bivariate_feat_2 = single_bivariate_feat;
single_bivariate_feat_2(incomplete_channels) = [];

for f = 1:5
    for r = [1:16,17:20]
        
        new_roi_2(new_roi_2==(2*r-1)) = 2*r;
        
        pred_uninvolved = all_feat_zscores2(find((class_ind==1).*(new_roi_2==2*r)),(f+5));
        pred_eiz = all_feat_zscores2(find((class_ind==2).*(new_roi_2==2*r)),(f+5));
        pred_soz = all_feat_zscores2(find((class_ind==3).*(new_roi_2==2*r)),(f+5));
        
        median(pred_uninvolved)
        median(pred_eiz)
        median(pred_soz)
        
        plot_matrix = padcat(pred_uninvolved, pred_eiz,pred_soz);
        
        p1_region(f,r) = ranksum(pred_uninvolved,pred_eiz)
        p2_region(f,r) = ranksum(pred_uninvolved,pred_soz)
        p3_region(f,r) = ranksum(pred_eiz,pred_soz)
        
        cohen_d2(f,r) = computeCohen_d(pred_uninvolved,pred_soz)
        cohen_d3(f,r) = computeCohen_d(pred_eiz,pred_soz)
        
    end
end

figure(1);clf;
subplot(3,1,1)
imagesc([p3_region])
xticks([1:20])
xticklabels(split_roi_name)
xtickangle(-45)
colorbar
subplot(3,1,2)
imagesc([p3_region]<0.05)
xticks([1:20])
xticklabels(split_roi_name)
xtickangle(-45)
subplot(3,1,3)
imagesc([cohen_d3])
xticks([1:20])
xticklabels(split_roi_name)
xtickangle(-45)
caxis([0,1])
colorbar

%% patient-level abnormality score res in/out good/poor outcome
good_rz = [];
good_out = [];

poor_rz = [];
poor_out = [];

for s = 107:166
    this_pt_outcome = late_outcome(s-106);
    this_pt_ind_rz = find([all_possible_labels(:,1)==s].*[all_possible_labels(:,5)==1]);
    this_pt_ind_out = find([all_possible_labels(:,1)==s].*[all_possible_labels(:,5)==0]);
    
    if this_pt_outcome==1
        good_rz = [good_rz; nanmedian(all_pred_analysis(this_pt_ind_rz))];
        good_out = [good_out; nanmedian(all_pred_analysis(this_pt_ind_out))];
    else
        poor_rz = [poor_rz; nanmedian(all_pred_analysis(this_pt_ind_rz))];
        poor_out = [poor_out; nanmedian(all_pred_analysis(this_pt_ind_out))];
    end
end

plot_matrix = padcat(good_rz,good_out,poor_rz,poor_out);
figure(2);clf;
UnivarScatter(plot_matrix)
