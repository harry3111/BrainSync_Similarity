clear
clc

setenv('LD_LIBRARY_PATH',[ '/usr/lib:/home/zhark2/abin:' getenv('LD_LIBRARY_PATH') ]) 
%setenv('LD_LIBRARY_PATH',[ '/usr/lib:/Users/salzwedelap/abin:' getenv('DYLD_LIBRARY_PATH') ]) 


%Path to AFNI directory.
%setenv('PATH', [getenv('PATH') ':/Users/salzwedelap/abin'])
setenv('PATH', [getenv('PATH') ':/home/zhark2/abin'])
% FSL related (link)
setenv( 'FSLDIR', '/usr/local/fsl' );
setenv( 'FSLOUTPUTTYPE', 'NIFTI_GZ');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);

%% ttest of matching MOC - other MOCs (use z-transformed corr!)  
 
%clust_thrs = {'1959.0', '3705.0'};  %% Infant's or MOC's; alpha=0.05(1-sided!!!!,NN=1), p=p_thr. FD='_FD0.5mm',threshold=210   %%!! 
%clust_thrs = {'1944.0', '2878.2'};  %% Infant's or MOC's; alpha=0.05(1-sided!!!!,NN=1), p=0.05. FD='_FD0.5mm',threshold=360  %%!!  
%clust_thrs = {'403.9', '532.7'};  %% Infant's or MOC's; alpha=0.05(1-sided!!!!,NN=1), p=0.01. FD='_FD0.5mm',threshold=360  %%!!  
clust_thrs = {'79.8', '103.0'};  %% Infant's or MOC's; alpha=0.05(1-sided!!!!,NN=1), p=0.001. FD='_FD0.5mm',threshold=360  %%!!  
FD = '_FD0.5mm'; %% 0.3mm:'' or 0.5mm:'_FD0.5mm' 
threshold = 360;  %%  200/210 or 360  
p_thr = '0.001';  %% 1-sided!!!!, RIGHT_TAIL 

clust_thr = num2str((str2double(clust_thrs{1}) + str2double(clust_thrs{2}))/2); % Infant & MOC
alldir = ['/media/zhark2/glab2/Haitao/New_DU_data_2/brainsync/BrainSync_Corr_result' FD '/tr' num2str(threshold) '_corr/']; 
%% transform Corr to Z 
datapath = [alldir];  
outputpath = [alldir 'ttest_matching-others/Z/'];  
if ~exist([outputpath 'matching_MOC'], 'dir')
    mkdir([outputpath 'matching_MOC'])
end
if ~exist([outputpath 'other_MOC'], 'dir')
    mkdir([outputpath 'other_MOC'])
end

subjects_Infant = importSubjIDs(['/media/zhark2/glab2/Haitao/New_DU_data_2/lists' FD '/Infant_list_gsr_tr' num2str(threshold) FD '_all_2.txt']);  %%!!  
num_Infant = length(subjects_Infant);
subjects_MOC = importSubjIDs(['/media/zhark2/glab2/Haitao/New_DU_data_2/lists' FD '/MOC_list_gsr_tr' num2str(threshold) FD '_all_2.txt']);  %%!!  
num_MOC = length(subjects_MOC);

for ss_Infant = 1:num_Infant
    subj_Infant = subjects_Infant{ss_Infant};
    infile_list = dir([datapath 'matching_MOC/' subj_Infant '_matchingMOC_*_brainsynced_CorrMap.nii.gz']); 
    if ~isempty(infile_list)
        % matching MOC
        infile = [datapath 'matching_MOC/' infile_list.name];
        outfile = [outputpath 'matching_MOC/' infile_list.name(1:end-7) '_Z.nii.gz'];
        if exist(infile,'file') && ~exist(outfile,'file')
            [~,~] = system(['3dcalc -a ' infile ' -expr ''atanh(a)'' -prefix ' outfile]);
        end
        % other MOCs
        datadir = [datapath 'other_MOC/' subj_Infant '_otherMOC'];
        outdir = [outputpath 'other_MOC/' subj_Infant '_otherMOC'];
        if ~exist(outdir, 'dir')
            mkdir(outdir);
        end
        for ss_MOC = 1:num_MOC
            subj_MOC = subjects_MOC{ss_MOC};
            infile = [datadir '/' subj_Infant '_otherMOC_' subj_MOC '_brainsynced_CorrMap.nii.gz'];
            outfile = [outdir '/' subj_Infant '_otherMOC_' subj_MOC '_brainsynced_CorrMap_Z.nii.gz'];
            if exist(infile,'file') && ~exist(outfile,'file')
                [~,~] = system(['3dcalc -a ' infile ' -expr ''atanh(a)'' -prefix ' outfile]);
            end
        end
        % make average Z map for other MOCs
        infiles_list = rdir([outdir '/' subj_Infant '_otherMOC_*_brainsynced_CorrMap_Z.nii.gz']);
        dset = [];
        for ff = 1:length(infiles_list)
            dset = [dset ' ' infiles_list(ff).name];
        end
        outfile = [outputpath 'other_MOC/' subj_Infant '_otherMOC_averaged_brainsynced_CorrMap_Z.nii.gz'];  %%!!  
        if ~exist(outfile, 'file')
            [~,~] = system(['3dMean -prefix ' outfile  dset]); 
        end
    end
end


%% ttest of matching MOC-others using z-transformed Z maps 
datapath = [alldir 'ttest_matching-others/Z/'];  
outputpath = [alldir 'ttest_matching-others/ttest_results/'];
if ~exist(outputpath, 'dir')
    mkdir(outputpath);
end

subjects_Infant = importSubjIDs(['/media/zhark2/glab2/Haitao/New_DU_data_2/lists' FD '/Infant_list_gsr_tr' num2str(threshold) FD '_all_2.txt']);  %%!!  
num_Infant = length(subjects_Infant);

for ss_Infant = 1:num_Infant
    subj_Infant = subjects_Infant{ss_Infant};
    infileA_list = rdir([datapath 'matching_MOC/' subj_Infant '_matchingMOC_*_brainsynced_CorrMap_Z.nii.gz']);
    if ~isempty(infileA_list)
        infileA = infileA_list.name;
        setB = [datapath 'other_MOC/' subj_Infant '_otherMOC/' subj_Infant '_otherMOC_*_brainsynced_CorrMap_Z.nii.gz'];
        outfile = [outputpath subj_Infant '_ttest_matchingMOC-otherMOCs.nii.gz'];
        if exist(infileA,'file') && ~exist(outfile,'file')
            system(['3dttest++ -singletonA ' infileA ' -setB ' setB ' -prefix ' outfile]); 
        end
        %% cluster-correction
        infile = [outputpath subj_Infant '_ttest_matchingMOC-otherMOCs.nii.gz'];
        outdir = [outputpath 'Cluster_corrected_p' p_thr '/'];
        if ~exist(outdir, 'dir')
            mkdir(outdir);
        end
        outfile = [outdir 'Cluster_corrected_p' p_thr '_' subj_Infant '_ttest_matchingMOC-otherMOCs.nii.gz'];
        outfile_1 = [outdir 'ClusterMask_Cluster_corrected_p' p_thr '_' subj_Infant '_ttest_matchingMOC-otherMOCs.nii.gz'];  %%!!  
        if exist(infile, 'file') && ~(exist(outfile, 'file') && exist(outfile_1, 'file'))
            system(['3dClusterize -nosum -1Dformat -inset ' infile ' -idat 1 -ithr 1 -NN 1 -clust_nvox ' clust_thr ' -1sided RIGHT_TAIL p=' p_thr ' -pref_dat ' outfile ' -pref_map ' outfile_1]); %% 
        end
    end
end

%% ttest of matching MOCs-averaged_other MOCs
datapath = [alldir 'ttest_matching-others/Z/'];  
outputpath = [alldir 'ttest_matching-others/Z/'];  
if ~exist(outputpath, 'dir')
    mkdir(outputpath);
end

setA = [datapath 'matching_MOC/*_matchingMOC_*_brainsynced_CorrMap_Z.nii.gz'];
setB = [datapath 'other_MOC/*_otherMOC_averaged_brainsynced_CorrMap_Z.nii.gz'];
outfile = [outputpath 'matchingMOC-otherMOC_averaged_ttest_paired.nii.gz'];
if ~exist(outfile, 'file')
    system(['3dttest++ -setA ' setA ' -setB ' setB ' -prefix ' outfile ' -paired']); 
end
%% cluster-correction
infile = [outputpath 'matchingMOC-otherMOC_averaged_ttest_paired.nii.gz'];
outfile = [outputpath 'Cluster_corrected_p' p_thr '_matchingMOC-otherMOC_averaged_ttest_paired.nii.gz'];
outfile_1 = [outputpath 'ClusterMask_Cluster_corrected_p' p_thr '_matchingMOC-otherMOC_averaged_ttest_paired.nii.gz'];  %%!!  
if exist(infile, 'file') && ~(exist(outfile, 'file') && exist(outfile_1, 'file'))
    system(['3dClusterize -nosum -1Dformat -inset ' infile ' -idat 1 -ithr 1 -NN 1 -clust_nvox ' clust_thr ' -1sided RIGHT_TAIL p=' p_thr ' -pref_dat ' outfile ' -pref_map ' outfile_1]); %% 
end

%% ttest of matching MOCs 
datapath = [alldir 'ttest_matching-others/Z/matching_MOC/'];
outputpath = [alldir 'ttest_matching-others/Z/matching_MOC/'];
if ~exist(outputpath, 'dir')
    mkdir(outputpath);
end

setA = [datapath '*_matchingMOC_*_brainsynced_CorrMap_Z.nii.gz'];
outfile = [outputpath 'matchingMOC_ttest.nii.gz'];
if ~exist(outfile, 'file')
    system(['3dttest++ -setA ' setA ' -prefix ' outfile]); 
end
%% cluster-correction
infile = [outputpath 'matchingMOC_ttest.nii.gz'];
outfile = [outputpath 'Cluster_corrected_p' p_thr '_matchingMOC_ttest.nii.gz'];
outfile_1 = [outputpath 'ClusterMask_Cluster_corrected_p' p_thr '_matchingMOC_ttest.nii.gz'];  %%!!  
if exist(infile, 'file') && ~(exist(outfile, 'file') && exist(outfile_1, 'file'))
    system(['3dClusterize -nosum -1Dformat -inset ' infile ' -idat 1 -ithr 1 -NN 1 -clust_nvox ' clust_thr ' -1sided RIGHT_TAIL p=' p_thr ' -pref_dat ' outfile ' -pref_map ' outfile_1]); %% 
end

%% ttest of averaged_other MOCs 
datapath = [alldir 'ttest_matching-others/Z/other_MOC/'];
outputpath = [alldir 'ttest_matching-others/Z/other_MOC/'];
if ~exist(outputpath, 'dir')
    mkdir(outputpath);
end

setA = [datapath '*_otherMOC_averaged_brainsynced_CorrMap_Z.nii.gz'];
outfile = [outputpath 'otherMOC_averaged_ttest.nii.gz'];
if ~exist(outfile, 'file')
    system(['3dttest++ -setA ' setA ' -prefix ' outfile]); 
end
%% cluster-correction
infile = [outputpath 'otherMOC_averaged_ttest.nii.gz'];
outfile = [outputpath 'Cluster_corrected_p' p_thr '_otherMOC_averaged_ttest.nii.gz'];
outfile_1 = [outputpath 'ClusterMask_Cluster_corrected_p' p_thr '_otherMOC_averaged_ttest.nii.gz'];  %%!!  
if exist(infile, 'file') && ~(exist(outfile, 'file') && exist(outfile_1, 'file'))
    system(['3dClusterize -nosum -1Dformat -inset ' infile ' -idat 1 -ithr 1 -NN 1 -clust_nvox ' clust_thr ' -1sided RIGHT_TAIL p=' p_thr ' -pref_dat ' outfile ' -pref_map ' outfile_1]); %% 
end

