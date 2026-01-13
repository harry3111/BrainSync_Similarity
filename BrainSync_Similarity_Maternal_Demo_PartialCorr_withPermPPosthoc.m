clear
clc

datapath = ['/media/zhark2/glab2/Haitao/New_DU_data_2/brainsync/BrainSync_Corr_result_FD0.5mm/tr360_corr/ttest_matching-others/Z/'];
outputpath = ['/media/zhark2/glab2/Haitao/New_DU_data_2/brainsync/'];
sim_str = 'BrainSync_Matched_FrontalCluster_Similarity';  %%!!  WholeBrain / FrontalCluster
sim_sstr = 'BrainSync-Matched-FrontalCluster-Similarity';  %%!!  WholeBrain / FrontalCluster
sim_file = [datapath 'matched_unmatched_Cluster_1_p0.001_S_ttest_results.mat'];  %%!!  wb / Cluster_1_p0.001
load(sim_file);
Sim = matched_cluster_S;  %%!!  wb / cluster
FD = '_FD0.5mm'; %% 0.3mm:'' or 0.5mm:'_FD0.5mm' 
threshold = 360;  %%  200/210 or 360  
group_1 = 'Infant';
group_2 = 'MOC';
subjects_matched_1 = importSubjIDs(['/media/zhark2/glab2/Haitao/New_DU_data_2/lists' FD '/' group_1 '_list_truepair_tr' num2str(threshold) '_updated.txt']);
subjects_matched_2 = importSubjIDs(['/media/zhark2/glab2/Haitao/New_DU_data_2/lists' FD '/' group_2 '_list_truepair_tr' num2str(threshold) '_updated.txt']);
outdir = [outputpath 'Sim_Demo_PartialCorr_Results_' sim_str '_Posthoc_CORRECTED/'];  %%!!    %%!!  
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

demo_dir = ['/media/zhark2/glab2/Haitao/New_DU_data_2/lists_FD0.5mm/'];
demo_file = [demo_dir 'Maternal_mood_table_truepair_20241121_CORRECTED.mat'];  %%!!  extracted demo    (updated!!)      %%!!  
load(demo_file);
demo_table = maternal_mood_table_truepair;  %%!!  
cov_strs = {'infant_sex_male', 'birthweight_grams', 'gestAge_birth_days', 'gestAge_scan_days', 'meanFDs_Infant_truepair', 'meanFDs_MOC_truepair', 'avgPreINR', 'MOC_educationYears'};  %%!!    %%!!  
cov_sstrs = {'infant-sex-male', 'birthweight-grams', 'gestAge-birth-days', 'gestAge-scan-days', 'meanFDs-Infant-truepair', 'meanFDs-MOC-truepair', 'avgPreINR', 'MOC-educationYears'};  %%!!    %%!!  
ind_valCovs = ones(size(demo_table, 1), 1, 'logical');
for cc = 1:(length(cov_strs))  %%!!    %%!!  
    cov_str = cov_strs{cc};
    cov_sstr = cov_sstrs{cc};
    eval(['Cov = demo_table.' cov_str ';']);
    ind_valCov = ~isnan(Cov); %% 
    ind_valCovs = ind_valCovs & ind_valCov;  %%!!    %%!!  
end
%load([demo_dir 'meanFDs_Infant_truepair.mat']); %% 
%meanFDs_Infant_truepair = meanFDs;
%load([demo_dir 'meanFDs_MOC_truepair.mat']); %% 
%meanFDs_MOC_truepair = meanFDs;
%demo_strs = {'CESD_preAVG', 'EPDS_preAVG', 'STATE_preAVG', 'CESD_tot_all_pv1', 'EPDS_tot_all_pv1', 'STATE_Tot_all_pv1', 'CESD_tot_all_pv2', 'EPDS_tot_all_pv2', 'STATE_Tot_all_pv2', 'CESD_tot_all_pv3', 'EPDS_tot_all_pv3', 'STATE_Tot_all_pv3', 'avgPreINR', 'MOC_educationYears', 'Sensitivity_Include', 'Structuring_Include', 'Nonintrusiveness_Include', 'Nonhostility_Include', 'total_problem_Include', 'total_competence_Include'};  %%!!  
%demo_sstrs = {'CESD-preAVG', 'EPDS-preAVG', 'STATE-preAVG', 'CESD-tot-all-pv1', 'EPDS-tot-all-pv1', 'STATE-Tot-all-pv1', 'CESD-tot-all-pv2', 'EPDS-tot-all-pv2', 'STATE-Tot-all-pv2', 'CESD-tot-all-pv3', 'EPDS-tot-all-pv3', 'STATE-Tot-all-pv3', 'avgPreINR', 'MOC-educationYears', 'Sensitivity-Include', 'Structuring-Include', 'Nonintrusiveness-Include', 'Nonhostility-Include', 'total-problem-Include', 'total-competence-Include'};  %%!!  
demo_strs = {'Sensitivity_Include', 'Structuring_Include', 'Nonintrusiveness_Include', 'Nonhostility_Include'};  %%!!    %%!!  
demo_sstrs = {'Sensitivity-Include', 'Structuring-Include', 'Nonintrusiveness-Include', 'Nonhostility-Include'};  %%!!    %%!!  
valPs = zeros(length(demo_strs), 1); %% 
valPermPs = zeros(length(demo_strs), 1); %% 
for dd = 1:length(demo_strs)
    demo_str = demo_strs{dd};
    demo_sstr = demo_sstrs{dd};
    eval(['Demo = demo_table.' demo_str ';']);
    ind_valDemo = ~isnan(Demo); %% 
    ind_valDemo = ind_valDemo & ind_valCovs;  %%!!    %%!!  
    valDemo = Demo(ind_valDemo);
    valSim = Sim(ind_valDemo);
    valCovs = [];
    for cc = 1:(length(cov_strs))  %%!!    %%!!  
        cov_str = cov_strs{cc};
        cov_sstr = cov_sstrs{cc};
        eval(['Cov = demo_table.' cov_str ';']);
        valCov = Cov(ind_valDemo);
        valCovs = [valCovs, valCov];  %%!!  
    end
    %Cov = meanFDs_Infant_truepair;
    %valCov = Cov(ind_valDemo);
    %valCovs = [valCovs, valCov];  %%!!  
    %Cov = meanFDs_MOC_truepair;
    %valCov = Cov(ind_valDemo);
    %valCovs = [valCovs, valCov];  %%!!  
    valsubjects_matched_1 = subjects_matched_1(ind_valDemo);
    valsubjects_matched_2 = subjects_matched_2(ind_valDemo);
    outfile1 = [outdir 'Sim_' sim_str '_Demo_' demo_str '_PartialCorr.mat'];
    outfile2 = [outdir 'Sim_' sim_str '_Demo_' demo_str '_PartialCorr.fig'];
    if ~(exist(outfile1, 'file') && exist(outfile2, 'file'))
        [valPartialCorr, valP] = partialcorr(valSim, valDemo, valCovs);  %%!!    %%!!  
        valPs(dd) = valP; %% 
        % Permutation P Values (valPermP)
        rand_seed = 1;  %% random seed  %%  
        perm_n = 1000;  %% number of permutations  %%  
        perm_valDemos = cell(perm_n, 1);  
        rng(rand_seed);
        for pp = 1:perm_n
            perm_valDemos{pp} = valDemo(randperm(length(valDemo))); 
        end
        perm_valPartialCorrs = zeros(perm_n, 1); 
        for pp = 1:perm_n
            perm_valPartialCorrs(pp) = partialcorr(valSim, perm_valDemos{pp}, valCovs); %% 
        end
        if valPartialCorr < 0
            smaller_cnt = sum(perm_valPartialCorrs < valPartialCorr); 
        elseif valPartialCorr >= 0
            smaller_cnt = sum(perm_valPartialCorrs > valPartialCorr);
        else
            smaller_cnt = NaN; 
        end
        valPermP = smaller_cnt / perm_n;  %%!!  
        valPermPs(dd) = valPermP; %% 
        save(outfile1, 'sim_str', 'sim_sstr', 'cov_strs', 'cov_sstrs', 'demo_str', 'demo_sstr', 'Sim', 'Demo', 'ind_valDemo', 'valSim', 'valDemo', 'valCovs', 'group_1', 'group_2', 'subjects_matched_1', 'subjects_matched_2', 'valsubjects_matched_1', 'valsubjects_matched_2', 'valPartialCorr', 'valP', 'perm_valDemos', 'perm_valPartialCorrs', 'valPermP', 'ind_valCovs');  %%!!  
        
        % scatter plots    
        x = valSim;
        y = valDemo;
        h_fig = figure(); 
        hold on;
        plt = plot(x, y, '.', 'MarkerSize', 20, 'Color', 'k');  %%  
        P = polyfit(x,y,1);
        x0 = min(x) ; x1 = max(x) ;
        xi = linspace(x0,x1) ;
        yi = P(1)*xi+P(2);
        l = plot(xi,yi,'k-') ;  %%  
        l.DisplayName = 'regression line';
        l.LineWidth = 2;
        xlabel('Similarity');
        ylabel('Demographic');
        title(['Sim ' sim_sstr ' - Demo ' demo_sstr]); %% 
        %xlim([ ]); %% 
        %ylim([-0.6 0.8]); %% 
        hold off;
        saveas(h_fig, outfile2);
        close(h_fig);  %  
        
    end
end
outfile = [outdir 'Ps_FDR_H_Results_' sim_str '.mat'];
if ~exist(outfile, 'file')
    alpha = 0.05; %% 
    valPs_FDR = mafdr(valPs, 'BHFDR', true); % few terms: use BHFDR 
    valPs_FDR_H = valPs_FDR <= alpha;
    valPermPs_FDR = mafdr(valPermPs, 'BHFDR', true); % few terms: use BHFDR 
    valPermPs_FDR_H = valPermPs_FDR <= alpha;
    valPs_H = valPs <= alpha;
    valPermPs_H = valPermPs <= alpha;
    save(outfile, 'alpha', 'valPs_FDR', 'valPs_FDR_H', 'valPermPs_FDR', 'valPermPs_FDR_H', 'valPs', 'valPs_H', 'valPermPs', 'valPermPs_H', 'cov_strs', 'demo_strs'); %% 
end

