% extra code snippets
%% this code can be used to visualize distributions of specific edges

res_conn = [];
for roi1 = 17
    for roi2 = 39
        
        % assemble all edges in hup patients
        for s = 1:max(all_pts)
            pt_roi_new = new_roi(all_pts==s);
            %pt_res = resect_ch_all(all_pts==s);
            res_roi = pt_roi_new;%(find(pt_res));
            
            if sum(res_roi==roi1)+sum(res_roi==roi2)
                
                % assemble the nodes together
                pt_adj = all_pt_adj{s};
                
                % extract resected but only if roi1 and roi2 affected
                adj_res = pt_adj(find(pt_roi_new==roi1),find(pt_roi_new==roi2));
                
                res_conn = [res_conn;adj_res(:)];
                
              
            end
        end



% conn 1-27
edge_hup_roi1_roi2 = raw_atlas_edge_hup(roi1,roi2,:);
edge_hup_roi1_roi2(isnan(edge_hup1_27)) = [];
edge_mni_roi1_roi2 = raw_atlas_edge_mni(roi1,roi2,:);
edge_mni_roi1_roi2(isnan(edge_mni1_27)) = [];

figure(2);clf;
hold on
histogram((squeeze(edge_mni_roi1_roi2)),10,'Normalization','probability')
histogram((squeeze(edge_hup_roi1_roi2)),10,'Normalization','probability')

hold off
        try all_conn_pval = ranksum(squeeze(edge_mni1_27),squeeze(edge_hup1_27))
        catch anyerror
            all_conn_pval = NaN;
        end
        
    end
end
%% entropy of specific edges
soz_ch_bin = [zeros(1772,1);soz_ch];

all_entropy = all_feat_zscores(:,6);
all_entropy_1 = all_entropy(find((new_roi==5).*soz_ch_bin));
all_entropy_2 = all_entropy(find((new_roi==26).*soz_ch_bin));
all_entropy_3 = all_entropy(find((new_roi==20).*soz_ch_bin));
all_entropy_4 = all_entropy(find((new_roi==17).*soz_ch_bin));
all_entropy_5 = all_entropy(find((new_roi==5).*(soz_ch_bin==0)));
all_entropy_6 = all_entropy(find((new_roi==26).*(soz_ch_bin==0)));
all_entropy_7 = all_entropy(find((new_roi==20).*(soz_ch_bin==0)));
all_entropy_8 = all_entropy(find((new_roi==17).*(soz_ch_bin==0)));

figure(1);clf;
subplot(2,2,1)
hold on
histogram(all_entropy_5,10,'Normalization','probability')
histogram(all_entropy_1,10,'Normalization','probability')
ranksum(all_entropy_5,all_entropy_1)
hold off
subplot(2,2,2)
hold on
histogram(all_entropy_6,10,'Normalization','probability')
histogram(all_entropy_2,10,'Normalization','probability')
ranksum(all_entropy_6,all_entropy_2)
hold off
subplot(2,2,3)
hold on
histogram(all_entropy_7,10,'Normalization','probability')
histogram(all_entropy_3,10,'Normalization','probability')
ranksum(all_entropy_7,all_entropy_3)
hold off
subplot(2,2,4)
hold on
histogram(all_entropy_8,10,'Normalization','probability')
histogram(all_entropy_4,10,'Normalization','probability')
ranksum(all_entropy_8,all_entropy_4)
hold off
%% new areas to test

figure(1);clf
for i = 1:21
    %roi_mni = intersect(all_normal_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_mni = intersect(normal_MNI_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    %roi_hup = intersect(EIZ_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_hup = intersect(normal_HUP_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    %roi_abn = intersect([soz_good_ch;soz_poor_ch],[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    
    % mean mni
    mean_mni = median(pxx_norm(:,roi_mni)');
    mean_hup = median(pxx_norm(:,roi_hup)');
    %mean_abn = median(pxx_norm(:,roi_abn)');
    
    num_ch_roi(i,1) = length(roi_mni);
    num_ch_roi(i,2) = length(roi_hup);
    %num_ch_roi(i,3) = length(roi_abn);

    prct_5_95_mni = prctile(pxx_norm(:,roi_mni)',[25,75]);
    prct_5_95_hup = prctile(pxx_norm(:,roi_hup)',[25,75]);
    %prct_5_95_abn = prctile(pxx_norm(:,roi_abn)',[25,75]);
    
    %region_percentile_mni(i).data = prct_5_95_mni;
    
    color1 = [0, 0.4470, 0.7410];
    color2 = [0.6350, 0.0780, 0.1840];
    color3 = [255, 248, 209]./255; 
    color4 = [103 55 155]/255;
    color6 = [78 172 91]/255;

    
    x = [0.5:0.5:80];
    
    subplot(3,7,i)
    hold on
    plot([4,4],[0,0.4],'k-')
    plot([8,8],[0,0.4],'k-')
    plot([13,13],[0,0.4],'k-')
    plot([30,30],[0,0.4],'k-')
    plot(x,mean_mni)
    patch([x fliplr(x)], [prct_5_95_mni(1,:) fliplr(prct_5_95_mni(2,:))], color1, 'FaceAlpha',0.5, 'EdgeColor','none')
    plot(x,mean_hup)
    patch([x fliplr(x)], [prct_5_95_hup(1,:) fliplr(prct_5_95_hup(2,:))], color6, 'FaceAlpha',0.5, 'EdgeColor','none')
    %plot(x,mean_abn)
%     try patch([x fliplr(x)], [prct_5_95_abn(1,:) fliplr(prct_5_95_abn(2,:))], color2, 'FaceAlpha',0.5, 'EdgeColor','none')
%     catch anyerror
%     end
%     plot(x,mean_mni)
%     patch([x fliplr(x)], [prct_5_95_mni(1,:) fliplr(prct_5_95_mni(2,:))], [0, .6, .77], 'FaceAlpha',0.5, 'EdgeColor','none')
%     %patch([x2 fliplr(x2)], [y2-CI2 fliplr(y2+CI2)], [.89, 0, .23], 'FaceAlpha',0.5, 'EdgeColor','none')
    %shadedErrorBar([0.5,2:80],mean_mni,sem_mni,'lineprops',{'markerfacecolor','red'})
    %shadedErrorBar([0.5,2:80],mean_hup,sem_hup,'lineprops',{'markerfacecolor','blue'})
%     try shadedErrorBar([0.5,2:80],mean_abn,sem_abn,'lineprops',{'markerfacecolor','green'})
%     catch anyerror
%     end
    this_roi_name = split(string(custom_atlas{2*i,1}),'_R');
    split_roi_name{i,1} = this_roi_name(1);
    title(sprintf('%s',this_roi_name(1)), 'Interpreter', 'none')
    xlim([0.5,80])
    ylim([0 0.2])
    
    % greek letters
    txt_d = '\delta';
    txt_t = '\theta';
    txt_a = '\alpha';
    txt_b = '\beta';
    txt_g = '\gamma';
    text(2,0.15,txt_d);
    text(5.5,0.15,txt_t);
    text(9.5,0.15,txt_a);
    text(18,0.15,txt_b);
    text(40,0.15,txt_g);
    
    set(gca, 'XScale','log', 'XTick',[0.5,4,8,13,30,80], 'XTickLabel',[0.5,4,8,13,30,80])
    ylabel('Spectral Density')
    xlabel('Frequency (Hz)')
    hold off
    
end