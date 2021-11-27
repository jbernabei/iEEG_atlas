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

all_entropy = all_pt_zscores(:,6);
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
%% check bivariate features in order to assess whether they distinguish
all_normal_ch = zeros(5203,1);
all_normal_ch(normal_HUP_ch) = 1;
all_normal_ch(normal_MNI_ch) = 1;
which_feat = 6;

% epileptogenic zone in the amygdala + hippocampus or not
MTL_channels = [new_roi==17]+[new_roi==18];
soz_ch_bin = [zeros(1772,1);soz_ch];
EIZ = zeros(5203,1);EIZ(EIZ_ch) = 1;
SOZ_MTL_channels = MTL_channels.*soz_ch_bin;%+MTL_channels.*EIZ;
EIZ_MTL_channels = MTL_channels.*EIZ;
norm_MTL_channels = MTL_channels.*all_normal_ch;

SOZ_MTL_feat = single_bivariate_feat(find(SOZ_MTL_channels));%,which_feat);
EIZ_MTL_feat = single_bivariate_feat(find(EIZ_MTL_channels));%,which_feat);
norm_MTL_feat = single_bivariate_feat(find(norm_MTL_channels));%,which_feat);

figure(1);clf;
hold on
scatter(ones(length(norm_MTL_feat),1),norm_MTL_feat,'jitter','on')
plot([0.75 1.25],[nanmedian(norm_MTL_feat) nanmedian(norm_MTL_feat)],'k-','LineWidth',2)
scatter(2*ones(length(EIZ_MTL_feat),1),EIZ_MTL_feat,'jitter','on')
plot([1.75 2.25],[nanmedian(EIZ_MTL_feat) nanmedian(EIZ_MTL_feat)],'k-','LineWidth',2)
scatter(3*ones(length(SOZ_MTL_feat),1),SOZ_MTL_feat,'jitter','on')
plot([2.75 3.25],[nanmedian(SOZ_MTL_feat) nanmedian(SOZ_MTL_feat)],'k-','LineWidth',2)
ylim([0,5])
hold off

p1 = ranksum(norm_MTL_feat,EIZ_MTL_feat)
p2 = ranksum(EIZ_MTL_feat,SOZ_MTL_feat)
p3 = ranksum(norm_MTL_feat,SOZ_MTL_feat)

%% get patient-level SOZ, RZ, spike indicators
all_normal_good = [];
all_normal_poor = [];
all_les_in_in = [];
all_nonles_in_in = [];
all_SOZ_in_in = [];
all_SOZ_in_out = [];
all_spike_in_in = [];
all_spike_in_out = [];

b = 0;
ptile = 75;
for s = 107:166
    b = b+1
    a = s;
    
    pt_roi = pt_loc{s};
    res = resected_ch([patient_no==b]); 
    spike = spike_24h([patient_no==b]);
    soz = soz_ch([patient_no==b]);
    
    roi_res = NaN*zeros(42,1);
    
%     for r = 1:42
%         this_roi = (r);
%         roi_res(r) = mean(res(pt_roi==this_roi));
%         roi_soz(r) = mean(soz(pt_roi==this_roi));
%         roi_spike(r) = mean(spike(pt_roi==this_roi));
%     end
%     
%     HUP_scores(a).res = ceil(roi_res);
%     HUP_scores(a).spike = (roi_spike>24);
%     HUP_scores(a).soz = ceil(roi_soz);
%     
%     res_roi = [HUP_scores(a).res==1];
%     soz_roi = ([HUP_scores(a).soz==1]);
%     spk_roi = ([HUP_scores(a).spike==1]);
%     
%     non_res_roi = ([HUP_scores(a).res==0]);
%     non_soz_roi = ([HUP_scores(a).soz==0]);
%     non_spk_roi = ([HUP_scores(a).spike==0]);
%     
    
    % spared edges
    all_normal = bivariate_native(feat).subj(a).data(find(~res),find(~res));
    HUP_normal(b,1) = nanmean(all_normal(:));
    
    all_normal_good = [all_normal_good;prctile((all_normal(:)),ptile)];
    
    % resected edges 
    %in_out_resect = bivariate_native(feat).subj(a).data(find(res),:);
    in_in_resect = bivariate_native(feat).subj(a).data(find(res),find(~res));
    %HUP_resect(b,1) = nanmedian(in_out_resect(:)); % (in-out)
    %HUP_resect(b,2) = nanmedian(in_in_resect(:)); % in-in
    %HUP_resect(b,3) = nanmedian(test_HUP_scores(:)); % all
    if lesion_field(b)
    all_les_in_in = [all_les_in_in;prctile(in_in_resect(:),ptile)];
    else
    all_nonles_in_in = [all_nonles_in_in;prctile(in_in_resect(:),ptile)];
    end
    %all_RZ_in_out = [all_RZ_in_out;prctile((in_out_resect(:)),ptile)];
    
    % soz edges (in-out)
    in_out_soz = bivariate_native(feat).subj(a).data(find(soz),find(~soz));
    in_in_soz = bivariate_native(feat).subj(a).data(find(soz),find(soz));
    %HUP_soz(b,1) = nanmedian(in_out_soz(:)); % (in-out)
    %HUP_soz(b,2) = nanmedian(in_in_soz(:)); % in-in
    
    all_SOZ_in_in = [all_SOZ_in_in;prctile(in_in_soz(:),ptile)];
    all_SOZ_in_out = [all_SOZ_in_out;prctile((in_out_soz(:)),ptile)];
    
    % spike edges 
    try in_out_spike = bivariate_native(feat).subj(a).data(find(spike),find(~spike));
        all_spike_in_out = [all_spike_in_out;prctile((in_out_spike(:)),ptile)];
    catch anyerror
    end
    try in_in_spike = bivariate_native(feat).subj(a).data(find(spike),find(spike));
        all_spike_in_in = [all_spike_in_in;prctile(in_in_spike(:),ptile)];
    catch anyerror
    end
%     HUP_spike(a,1) = nanmedian(in_out_spike(:)); % (in-out)
%     HUP_spike(a,2) = nanmedian(in_in_spike(:)); % in-in
%     
    
     
%     
end

%
figure(1);clf;
hold on

scatter(ones(length(all_normal_good),1),all_normal_good,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[nanmedian(all_normal_good) nanmedian(all_normal_good)],'k-','LineWidth',2)
scatter(3*ones(length(all_SOZ_in_in),1),all_SOZ_in_in,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([2.75 3.25],[nanmedian(all_SOZ_in_in) nanmedian(all_SOZ_in_in)],'k-','LineWidth',2)
scatter(4*ones(length(all_SOZ_in_out),1),all_SOZ_in_out,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([3.75 4.25],[nanmedian(all_SOZ_in_out) nanmedian(all_SOZ_in_out)],'k-','LineWidth',2)
scatter(5*ones(length(all_spike_in_in),1),all_spike_in_in,'MarkerEdgeColor',color4,'MarkerFaceColor',color4,'jitter','on')
plot([4.75 5.25],[nanmedian(all_spike_in_in) nanmedian(all_spike_in_in)],'k-','LineWidth',2)
scatter(6*ones(length(all_spike_in_out),1),all_spike_in_out,'MarkerEdgeColor',color4,'MarkerFaceColor',color4,'jitter','on')
plot([5.75 6.25],[nanmedian(all_spike_in_out) nanmedian(all_spike_in_out)],'k-','LineWidth',2)
xlim([0.5, 6.5])
ylim([0 6])
hold off

ranksum(all_normal_good,all_SOZ_in_in)
ranksum(all_normal_good,all_SOZ_in_out)
ranksum(all_normal_good,all_spike_in_in)
ranksum(all_normal_good,all_spike_in_out)

ranksum(all_SOZ_in_in,all_SOZ_in_out)
ranksum(all_SOZ_in_in,all_spike_in_in)
ranksum(all_spike_in_in,all_spike_in_out)
ranksum(all_SOZ_in_out,all_spike_in_out)
%% boxplots baby

num_good = sum([HUP_outcome==1]);
num_poor = sum([HUP_outcome>1]);

good_out=HUP_normal(HUP_outcome==1);
poor_out=HUP_normal(HUP_outcome>1);

good_in_out_soz=HUP_soz([HUP_outcome==1],1);
good_in_in_soz=HUP_soz([HUP_outcome==1],2);
poor_in_out_soz=HUP_soz([HUP_outcome>1],1);
poor_in_in_soz=HUP_soz([HUP_outcome>1],2);

good_in_out_spike=HUP_spike([HUP_outcome==1],1);
good_in_in_spike=HUP_spike([HUP_outcome==1],2);
poor_in_out_spike=HUP_spike([HUP_outcome>1],1);
poor_in_in_spike=HUP_spike([HUP_outcome>1],2);

all_plot_data = NaN*zeros(57,6);
all_plot_data(1:num_good,1) = good_in_in_soz;
all_plot_data(1:num_poor,2) = poor_in_in_soz;
all_plot_data(1:num_good,3) = good_in_out_soz;
all_plot_data(1:num_poor,4) = poor_in_out_soz;
all_plot_data(1:num_good,5) = good_out;
all_plot_data(1:num_poor,6) = poor_out;

ranksum(good_in_in_soz,poor_in_in_soz)
ranksum(good_in_out_soz,poor_in_out_soz)
ranksum(good_out,poor_out)

signrank(good_in_in_soz,good_in_out_soz)
signrank(good_in_out_soz,good_out)
signrank(good_in_in_soz,good_out)

signrank(poor_in_in_soz,poor_in_out_soz)
signrank(poor_in_out_soz,poor_out)
signrank(poor_in_in_soz,poor_out)

ranksum(good_in_in_spike,poor_in_in_spike)
ranksum(good_in_out_spike,poor_in_out_spike)
ranksum(good_out,poor_out)

signrank(good_in_in_spike,good_in_out_spike)
signrank(good_in_out_spike,good_out)
signrank(good_in_in_spike,good_out)

signrank(poor_in_in_spike,poor_in_out_spike)
signrank(poor_in_out_spike,poor_out)
signrank(poor_in_in_spike,poor_out)

figure(1);clf;
boxplot(all_plot_data)
ylim([-2 2])

%% lesional vs non-lesional
r1 = HUP_soz((lesion_field==0),2);
r2 = HUP_soz(lesion_field,2);

r3 = HUP_normal((lesion_field==0));
r4 = HUP_normal(lesion_field);

ranksum(r1,r2)

figure(1);clf;
hold on
plot(ones(length(r1),1),r1,'ko')
plot(2*ones(length(r2),1),r2,'ko')
xlim([0.5 2.5])
hold off


r1 = all_mean_wentropy([normal_HUP_ch;normal_MNI_ch]);
r2 = all_mean_wentropy(poor_out_ch);
r3 = all_mean_wentropy(EIZ_good_ch);
r4 = all_mean_wentropy(EIZ_poor_ch);
r5 = all_mean_wentropy(soz_good_ch);
r6 = all_mean_wentropy(soz_poor_ch);

ranksum(all_mean_wentropy(normal_HUP_ch),all_mean_wentropy(normal_MNI_ch))

p1 = ranksum(r1,r3)
p2 = ranksum(r3,r5)
p3 = ranksum(r1,r5)
p4 = ranksum(r2,r4)
p5 = ranksum(r4,r6)
p6 = ranksum(r2,r6)
p7 = ranksum(r1,r2)
p8 = ranksum(r3,r4)
p9 = ranksum(r5,r6)

figure(1);clf;
hold on
scatter(ones(length(r1),1),r1,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([0.75 1.25],[median(r1) median(r1)],'k-','LineWidth',2)
scatter(2*ones(length(r2),1),r2,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([1.75 2.25],[median(r2) median(r2)],'k-','LineWidth',2)
scatter(3*ones(length(r3),1),r3,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([2.75 3.25],[median(r3) median(r3)],'k-','LineWidth',2)
scatter(4*ones(length(r4),1),r4,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([3.75 4.25],[median(r4) median(r4)],'k-','LineWidth',2)
scatter(5*ones(length(r5),1),r5,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
plot([4.75 5.25],[median(r5) median(r5)],'k-','LineWidth',2)
scatter(6*ones(length(r6),1),r6,'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'jitter','on')
plot([5.75 6.25],[median(r6) median(r6)],'k-','LineWidth',2)
xlim([0.5, 6.5])
ylim([2,22])
ylabel('Wavelet entropy')
xticks([1:6])
xticklabels({'Engel 1, non-SOZ, non-EIZ','Engel 2+, non-SOZ, non-EIZ', 'Engel 1 EIZ','Engel 2+ EIZ','Engel 1 SOZ','Engel 2+ SOZ'})
hold off

%% wavelet entropy validation across all regions
for i = 1:60
    for j = 1:size(all_wake_data,2)
    start_inds = (i-1)*200+1;
    end_inds = 200*i;
    all_wentropy(i,j) = (wentropy(all_wake_data((start_inds:end_inds),j),'shannon'));
    end
end

all_mean_wentropy = log(-1*median(all_wentropy));


figure(1);clf;
for i = 1:20
    % extract indices that correspond to a given ROI for HUP and MNI data
    % (combining data from L and R hemispheres)
    roi_mni = intersect(normal_MNI_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_hup = intersect(normal_HUP_ch,[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_EIZ = intersect([EIZ_good_ch;EIZ_poor_ch],[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    roi_soz = intersect([soz_good_ch;soz_poor_ch],[find(new_roi==(2*(i-1)+1));find(new_roi==2*i)]);
    
    entropy_normal = all_mean_wentropy([roi_mni;roi_hup]);
    entropy_EIZ = all_mean_wentropy(roi_EIZ);
    entropy_soz = all_mean_wentropy(roi_soz);
    
    pval_EIZ(i,1) = ranksum(entropy_normal,entropy_EIZ);
    pval_soz(i,1) = ranksum(entropy_normal,entropy_soz);
    pval_EIZ_vs_soz(i,1) = ranksum(entropy_EIZ,entropy_soz);

    % median power spectral density across channels
    %entropy_mni = median(pxx_norm(:,roi_mni),2);
    %median_hup = median(pxx_norm(:,roi_hup),2);
    plot_cell{1,1} = entropy_normal;
    plot_cell{1,2} = entropy_EIZ;
    plot_cell{1,3} = entropy_soz;
    
    subplot(4,5,i)
    hold on
    violin(plot_cell,'xlabel',{'','',''},'facecolor',[color1;color7;color2],'mc',[],'medc','k');%,'edgecolor','k');
    legend('off')
    txt_sig1 = '+';
    txt_sig2 = '*';
    txt_sig3 = '#';
    ylim([7 23])
    
    if pval_EIZ(i,1) < (0.05./20)
        plot([1,2],[20,20],'k-','LineWidth',2)
        text(1.45,20.25,txt_sig2,'FontSize', 20);
    end
    if pval_soz(i,1) < (0.05./20)
        plot([1,3],[21.75,21.75],'k-','LineWidth',2)
        text(2,22,txt_sig2,'FontSize', 20);
    end
    if pval_EIZ_vs_soz(i,1) < (0.05./20)
        plot([2,3],[20,20],'k-','LineWidth',2)
        text(2.5,20.25,txt_sig2,'FontSize', 20);
    end
    
    ylabel('Log -Entropy')
    this_roi_name = split(string(custom_atlas{2*i,1}),'_R'); % extract roi name without laterality
    split_roi_name{i,1} = this_roi_name(1);
    title(sprintf('%s',this_roi_name(1)), 'Interpreter', 'none') % use that as title
    
    hold off
   
end

%% 
single_bivariate_feat_2 = single_bivariate_feat;
new_roi_2 = new_roi;
new_roi_2(incomplete_channels) = [];

new_roi_2(new_roi_2==17) = 18;

for ch = 1:4717
    if multi_class(ch)==1
        grouping{ch,1} = 'uninvolved';
        class_ind(ch,1) = 1;
    elseif multi_class(ch)==2
        grouping{ch,1} = 'EIZ';
        class_ind(ch,1) = 2;
    else
        grouping{ch,1} = 'SOZ';
        class_ind(ch,1) = 3;
    end
end

univariate_feats = all_feat_zscores(:,1:5);
univariate_feats(incomplete_channels,:) = [];

single_univariate_feat = max(abs(univariate_feats'))';

single_univariate_feat(isinf(single_univariate_feat)) = nanmedian(single_univariate_feat);
single_bivariate_feat_2(isinf(single_bivariate_feat_2)) = nanmedian(single_bivariate_feat_2);

single_univariate_feat = single_univariate_feat(new_roi_2==18);
single_bivariate_feat_2 = single_bivariate_feat_2(new_roi_2==18);
grouping = grouping(new_roi_2==18);
class_ind = class_ind(new_roi_2==18);


figure(1);clf;
h = scatterhist(single_univariate_feat,single_bivariate_feat_2,'Group',grouping,'Marker','.','MarkerSize',12)

clr = get(h(1),'colororder');
boxplot(h(2),single_univariate_feat,grouping,'orientation','horizontal',...
     'label',{'','',''},'color',clr);
      xlim(h(2),[0,30])
boxplot(h(3),single_bivariate_feat_2,grouping,'orientation','horizontal',...
     'label', {'','',''},'color',clr);
        xlim(h(3),[0,30])
 set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto'); 
xlim([0 10])
ylim([0 10])

d1 = computeCohen_d(single_univariate_feat(find(class_ind==2)),single_univariate_feat(find(class_ind==1)))
p1 = ranksum(single_univariate_feat(find(class_ind==2)),single_univariate_feat(find(class_ind==1)))
d2 = computeCohen_d(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==2)))
p2 = ranksum(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==2)))
d3 = computeCohen_d(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==1)))
p3 = ranksum(single_univariate_feat(find(class_ind==3)),single_univariate_feat(find(class_ind==1)))

d4 = computeCohen_d(single_bivariate_feat_2(find(class_ind==2)),single_bivariate_feat_2(find(class_ind==1)))
p4 = ranksum(single_bivariate_feat_2(find(class_ind==2)),single_bivariate_feat_2(find(class_ind==1)))
d5 = computeCohen_d(single_bivariate_feat_2(find(class_ind==3)),single_bivariate_feat_2(find(class_ind==2)))
p5 = ranksum(single_bivariate_feat_2(find(class_ind==3)),single_bivariate_feat_2(find(class_ind==2)))
d6 = computeCohen_d(single_bivariate_feat_2(find(class_ind==3)),single_bivariate_feat_2(find(class_ind==1)))
p6 = ranksum(single_bivariate_feat_2(find(class_ind==3)),single_bivariate_feat_2(find(class_ind==1)))

[d1 d2 d3 d4 d5 d6]
[p1 p2 p3 p4 p5 p6]
%% do node abormality (max univariate+bivariate in plot) for following:
% 1. within/outside soz for lesional / non-lesional

all_lesion = zeros(5203,1);
for i = 1:166
    if i>106
    all_lesion(all_pts==i) = lesion_field(i-106);
    end
end
 
within_soz_lesional = find(all_lesion.*soz_ch_all);
within_soz_nonlesional = find([all_lesion==0].*soz_ch_all);
outside_soz_lesional = find(all_lesion.*[soz_ch_all==0]);
outside_soz_nonlesional = find([all_lesion==0].*[soz_ch_all==0]);

for ch = 1:5203
    if ismember(ch,within_soz_lesional)
        lesional_grouping{ch,1} = 'within_soz_lesional';
    elseif ismember(ch,within_soz_nonlesional)
        lesional_grouping{ch,1} = 'within_soz_nonlesional';
    elseif ismember(ch,outside_soz_lesional)
        lesional_grouping{ch,1} = 'outside_soz_lesional';
    elseif ismember(ch,outside_soz_nonlesional)
        lesional_grouping{ch,1} = 'outside_soz_nonlesional';
    end
end

univariate_feats = all_feat_zscores(:,1:10);
single_univariate_feat = prctile(abs(univariate_feats'),100)';
single_bivariate_feat = prctile(abs_bivariate_feats',100)';

uni_within_soz_lesional = single_univariate_feat(within_soz_lesional);
uni_within_soz_nonlesional = single_univariate_feat(within_soz_nonlesional);
uni_outside_soz_lesional = single_univariate_feat(outside_soz_lesional);
uni_outside_soz_nonlesional = single_univariate_feat(outside_soz_nonlesional);
bi_within_soz_lesional = single_bivariate_feat(within_soz_lesional);
bi_within_soz_nonlesional = single_bivariate_feat(within_soz_nonlesional);
bi_outside_soz_lesional = single_bivariate_feat(outside_soz_lesional);
bi_outside_soz_nonlesional = single_bivariate_feat(outside_soz_nonlesional);

p1 = ranksum(uni_within_soz_lesional,uni_within_soz_nonlesional);
p2 = ranksum(uni_within_soz_lesional,uni_outside_soz_lesional);
p3 = ranksum(uni_within_soz_nonlesional,uni_outside_soz_nonlesional);
p4 = ranksum(uni_outside_soz_lesional,uni_outside_soz_nonlesional);
p5 = ranksum(bi_within_soz_lesional,bi_within_soz_nonlesional);
p6 = ranksum(bi_within_soz_lesional,bi_outside_soz_lesional);
p7 = ranksum(bi_within_soz_nonlesional,bi_outside_soz_nonlesional);
p8 = ranksum(bi_outside_soz_lesional,bi_outside_soz_nonlesional);

single_univariate_feat(incomplete_channels) = [];
single_bivariate_feat(incomplete_channels) = [];
lesional_grouping(incomplete_channels) = [];

figure(2);clf;
h = scatterhist(single_univariate_feat,single_bivariate_feat,'Group',lesional_grouping,'Marker','.','MarkerSize',12)

clr = get(h(1),'colororder');
boxplot(h(2),single_univariate_feat,lesional_grouping,'orientation','horizontal',...
     'label',{'','','',''},'color',clr);
      xlim(h(2),[0,30])
boxplot(h(3),single_bivariate_feat,lesional_grouping,'orientation','horizontal',...
     'label', {'','','',''},'color',clr);
        xlim(h(3),[0,30])
 set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto'); 
xlim([0 10])
ylim([0 10])

[p1 p2 p3 p4 p5 p6 p7 p8]
%% 2. within/outside ablation zone for good / poor outcome
all_outcome = ones(5203,1);
all_abl = NaN(5203,1);
for i = 1:166
    if i>106
    all_outcome(all_pts==i) = late_outcome(i-106);
        if therapy_field(i-106)
            all_abl(all_pts==i) = resect_ch_all(all_pts==i);
        end
    end
end

within_abl_good = find([all_outcome==1].*[all_abl==1]);
outside_abl_good = find([all_outcome==1].*[all_abl==0]);
within_abl_poor = find([all_outcome>1].*[all_abl==1]);
outside_abl_poor = find([all_outcome>1].*[all_abl==0]);

incomplete_channels2 = incomplete_channels;

abl_grouping = cell(5203,1);

for ch = 1:5203
    if ismember(ch,within_abl_good)
        abl_grouping{ch,1} = 'within_abl_good';
    elseif ismember(ch,outside_abl_good)
        abl_grouping{ch,1} = 'outside_abl_good';
    elseif ismember(ch,within_abl_poor)
        abl_grouping{ch,1} = 'within_abl_poor';
    elseif ismember(ch,outside_abl_poor)
        abl_grouping{ch,1} = 'outside_abl_poor';
    else
        incomplete_channels2 = [incomplete_channels2, ch];
    end
end

incomplete_channels2 = unique(incomplete_channels2);

univariate_feats = all_feat_zscores(:,1:10);
single_univariate_feat = prctile(abs(univariate_feats'),100)';
single_bivariate_feat = prctile(abs_bivariate_feats',100)';

uni_within_abl_good = single_univariate_feat(within_abl_good);
uni_outside_abl_good = single_univariate_feat(outside_abl_good);
uni_within_abl_poor = single_univariate_feat(within_abl_poor);
uni_outside_abl_poor = single_univariate_feat(outside_abl_poor);
bi_within_abl_good = single_bivariate_feat(within_abl_good);
bi_outside_abl_good = single_bivariate_feat(outside_abl_good);
bi_within_abl_poor = single_bivariate_feat(within_abl_poor);
bi_outside_abl_poor = single_bivariate_feat(outside_abl_poor);

p1 = ranksum(uni_within_abl_good,uni_outside_abl_good);
p2 = ranksum(uni_within_abl_poor,uni_outside_abl_poor);
p3 = ranksum(uni_within_abl_good,uni_within_abl_poor);
p4 = ranksum(uni_outside_abl_good,uni_outside_abl_poor);
p5 = ranksum(bi_within_abl_good,bi_outside_abl_good);
p6 = ranksum(bi_within_abl_poor,bi_outside_abl_poor);
p7 = ranksum(bi_within_abl_good,bi_within_abl_poor);
p8 = ranksum(bi_outside_abl_good,bi_outside_abl_poor);

single_univariate_feat(incomplete_channels2) = [];
single_bivariate_feat(incomplete_channels2) = [];
abl_grouping(incomplete_channels2) = [];

figure(3);clf;
h = scatterhist(single_univariate_feat,single_bivariate_feat,'Group',abl_grouping,'Marker','.','MarkerSize',12)

clr = get(h(1),'colororder');
boxplot(h(2),single_univariate_feat,abl_grouping,'orientation','horizontal',...
     'label',{'','','',''},'color',clr);
      xlim(h(2),[0,30])
boxplot(h(3),single_bivariate_feat,abl_grouping,'orientation','horizontal',...
     'label', {'','','',''},'color',clr);
        xlim(h(3),[0,30])
 set(h(2:3),'XTickLabel','');
view(h(3),[270,90]);  % Rotate the Y plot
axis(h(1),'auto'); 
xlim([0 10])
ylim([0 10])

[p1 p2 p3 p4 p5 p6 p7 p8]
