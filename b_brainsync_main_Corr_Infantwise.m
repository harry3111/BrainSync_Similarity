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

%% brainsync main for voxel-wise Corr maps(Infantwise)

FD = '_FD0.5mm'; %% 0.3mm:'' or 0.5mm:'_FD0.5mm' 
threshold = 360;  %%  200/210 or 360  

addpath('/media/zhark2/glab1/Haitao/Toolbox_codes/NIfTI_20140122/');
addpath('/media/zhark2/glab2/Haitao/New_DU_data_2/brainsync/BrainSync_Matlab/');

datapath = ['/media/zhark2/glab2/Haitao/New_DU_data_2/brainsync/BrainSync_Voxel' FD '/tr' num2str(threshold) '/'];   
outputpath = ['/media/zhark2/glab2/Haitao/New_DU_data_2/brainsync/BrainSync_Corr_result' FD '/tr' num2str(threshold) '_corr/'];   
if ~exist(outputpath, 'dir')
    mkdir(outputpath); 
end 
matching_form = 'RISE_motherinfant_numbers_updated.xlsx';  %%!!  
matching_info = readcell(matching_form);
matching_MOCs = matching_info(:,1);
matching_Infants = matching_info(:,2);
if ~exist([outputpath 'matching_MOC'], 'dir')
    mkdir([outputpath 'matching_MOC'])
end
if ~exist([outputpath 'other_MOC'], 'dir')
    mkdir([outputpath 'other_MOC'])
end

load('/media/zhark2/glab2/Haitao/New_DU_data_2/brainsync/Generate_ConnV/IndL_connV_MNI_2mm.mat');
IndL = IndL_aal;  %% cross_MNI_aal_2mm_mask.nii 

nii = load_nii('/media/zhark2/glab1/Haitao/DUdata/seed_corr_maps/MOC/auditoryNetwork/M80300030_APPA_merged_trc_gsr_auditoryNetwork_Z_000_INDIV/WB_Z_ROI_001.nii.gz'); 
Atlas = nii.img; %%

subjects_Infant = importSubjIDs(['/media/zhark2/glab2/Haitao/New_DU_data_2/lists' FD '/Infant_list_gsr_tr' num2str(threshold) FD '_all_2.txt']);  %%!!  
num_Infant = length(subjects_Infant);
subjects_MOC = importSubjIDs(['/media/zhark2/glab2/Haitao/New_DU_data_2/lists' FD '/MOC_list_gsr_tr' num2str(threshold) FD '_all_2.txt']);  %%!!  
num_MOC = length(subjects_MOC);

if ~exist(['lists' FD], 'dir')
    mkdir(['lists' FD]); 
end
fileID_Infant = fopen(['lists' FD '/Infant_list_truepair_tr' num2str(threshold) '_updated_check.txt'], 'w');  %%!!     
fileID_MOC = fopen(['lists' FD '/MOC_list_truepair_tr' num2str(threshold) '_updated_check.txt'], 'w');  %%!!     
for ss_Infant = 1:num_Infant
    subj_Infant = subjects_Infant{ss_Infant};
    matching_ind = find(ismember(matching_Infants,subj_Infant));
    if ~isempty(matching_ind)
        matching_MOC = matching_MOCs{matching_ind};
        infile_Infant = [datapath 'Infant/' subj_Infant '_Infant_MNI_aal_2mm_crossmasked.mat'];
        infile_mMOC = [datapath 'MOC/' matching_MOC '_MOC_MNI_aal_2mm_crossmasked.mat'];
        if exist(infile_Infant, 'file') && exist(infile_mMOC, 'file')
            fprintf(fileID_Infant, '%s\n', subj_Infant);
            fprintf(fileID_MOC, '%s\n', matching_MOC);
            % matching MOC 
            outfile_Corr_mMOC = [outputpath 'matching_MOC/' subj_Infant '_matchingMOC_' matching_MOC '_brainsynced_CorrMap.nii.gz'];
            if ~exist(outfile_Corr_mMOC, 'file')
                clear('subV_nor');
                load(infile_Infant);
                sub_i = subV_nor;
                clear('subV_nor');
                load(infile_mMOC);
                sub_j = subV_nor;

                ref = sub_j;
                [subSynced, O] = brainSync(ref, sub_i);
                [t, n] = size(subSynced);
                CorrMap_i_j = zeros(n, 1);
                for i = 1:n
                    CorrMap_i_j(i) = corr(subSynced(:, i), ref(:, i));
                end

                W=zeros(size(Atlas,1)*size(Atlas,2)*size(Atlas,3),1);
                W(IndL,:)=CorrMap_i_j;
                W=reshape(W,[size(Atlas,1),size(Atlas,2),size(Atlas,3)]);
                nii.img=W;
                
                save_nii(nii, outfile_Corr_mMOC);
            end
            
            % other MOCs 
            outdir_Corr_oMOC = [outputpath 'other_MOC/' subj_Infant '_otherMOC'];
            if ~exist(outdir_Corr_oMOC, 'dir')
                mkdir(outdir_Corr_oMOC);
            end
            for ss_MOC = 1:num_MOC
                subj_MOC = subjects_MOC{ss_MOC};
                if ~strcmp(subj_MOC, matching_MOC)
                    other_MOC = subj_MOC;
                    infile_oMOC = [datapath 'MOC/' other_MOC '_MOC_MNI_aal_2mm_crossmasked.mat'];
                    outfile_Corr_oMOC = [outdir_Corr_oMOC '/' subj_Infant '_otherMOC_' other_MOC '_brainsynced_CorrMap.nii.gz'];
                    if exist(infile_oMOC, 'file') && ~exist(outfile_Corr_oMOC, 'file')
                        clear('subV_nor');
                        load(infile_Infant);
                        sub_i = subV_nor;
                        clear('subV_nor');
                        load(infile_oMOC);
                        sub_j = subV_nor;

                        ref = sub_j;
                        [subSynced, O] = brainSync(ref, sub_i);
                        [t, n] = size(subSynced);
                        CorrMap_i_j = zeros(n, 1);
                        for i = 1:n
                            CorrMap_i_j(i) = corr(subSynced(:, i), ref(:, i));
                        end

                        W=zeros(size(Atlas,1)*size(Atlas,2)*size(Atlas,3),1);
                        W(IndL,:)=CorrMap_i_j;
                        W=reshape(W,[size(Atlas,1),size(Atlas,2),size(Atlas,3)]);
                        nii.img=W;
                        
                        save_nii(nii, outfile_Corr_oMOC);
                    end
                end
            end
            % make average map for other MOCs
            files_Corr_oMOC = rdir([outdir_Corr_oMOC '/' subj_Infant '_otherMOC_*_brainsynced_CorrMap.nii.gz']);
            dset = [];
            for ff = 1:length(files_Corr_oMOC)
                dset = [dset ' ' files_Corr_oMOC(ff).name];
            end
            outfile_avg_Corr_oMOC = [outputpath 'other_MOC/' subj_Infant '_otherMOC_averaged_brainsynced_CorrMap.nii.gz'];  %%!!  move these old files!!  
            if ~exist(outfile_avg_Corr_oMOC, 'file')
                [~,~] = system(['3dMean -prefix ' outfile_avg_Corr_oMOC  dset]); 
            end
        end
    end
end
fclose(fileID_Infant);
fclose(fileID_MOC);

