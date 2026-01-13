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

%% pre-steps before brainsync analysis 

FD = '_FD0.5mm'; %% 0.3mm:'' or 0.5mm:'_FD0.5mm' 
threshold = 360;  %%  200/210 or 360  

addpath('/media/zhark2/glab1/Haitao/Toolbox_codes/NIfTI_20140122/');
addpath('/media/zhark2/glab2/Haitao/New_DU_data_2/brainsync/BrainSync_Matlab/');

datapath = ['/media/zhark2/glab2/Haitao/New_DU_data_2/PreProResults' FD '/tr' num2str(threshold) '/']; 
outputpath = ['/media/zhark2/glab2/Haitao/New_DU_data_2/brainsync/BrainSync_Voxel' FD '/tr' num2str(threshold) '/'];
if ~exist(outputpath, 'dir')
    mkdir(outputpath); 
end
groups = {'Infant', 'MOC'};
n_g = length(groups);

load('/media/zhark2/glab2/Haitao/New_DU_data_2/brainsync/Generate_ConnV/IndL_connV_MNI_2mm.mat');
IndL = IndL_aal;  %% cross_MNI_aal_2mm_mask.nii 
connV = connV_aal;  %% 32mm 
h_var = 0.72;  %% changable

for gg = 1:n_g 
    group = groups{gg}; 
    subjects = importSubjIDs(['/media/zhark2/glab2/Haitao/New_DU_data_2/lists' FD '/' group '_list_gsr_tr' num2str(threshold) FD '_all_2.txt']); 
    num_subj = length(subjects);
    if ~exist([outputpath group],'dir')
        mkdir([outputpath group]);
    end
    
    for ss = 1:num_subj
        subj = subjects{ss};
        infile = [datapath group '_trc_gsr_MNI_masked/cf05_' subj '_rfMRI_REST_APPA_merged_SC_RegBP_sc_WS_trc_gsr_MNI_masked.nii.gz'];
        outfile = [outputpath group '/' subj '_' group '_MNI_aal_2mm_crossmasked.mat']; 
        if exist(infile,'file') && ~exist(outfile,'file')       
            nii = load_nii(infile);
            SubTC = nii.img;
            subm = reshape(SubTC, [size(SubTC,1)*size(SubTC,2)*size(SubTC,3), size(SubTC,4)]);
            subm = subm(IndL, :);
            subV_orig = subm';
            subm(isnan(subm(:))) = 0;
            [subm_tNLM, s, mTC, sTC, ~] = tNLM_SpatioTemporalData(subm, connV, h_var);  %space*time
            [subm, mean_vector, norm_vector] = normalizeData(subm_tNLM');  %time*space
            subV_nor = subm;
            
            save(outfile, 'subV_orig', 'subV_nor'); 
            
        end
    end
end

    
