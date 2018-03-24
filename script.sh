#KaWa July 2017 - ASL pipeline
#Copyright (c) 2017 Karolina Wartolowska

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:  The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#!/bin/sh

SUBJECTS_DIR=/vols/Data
OUTPUT_DIR=/vols/Output
MASK_DIR=/vols/Data/masks


for visit in bl hc;  do

for subject in 001 002 003 004 005; do

for cond in up down; do
 
### copy the "structural to standard space" transforms from vFREY.feat resuts
imcp $SUBJECTS_DIR/${visit}${subject}/strucfiles.anat/T1_biascorr_brain.nii.gz $SUBJECTS_DIR/${visit}${subject}/ASL_${cond}/struc_brain.nii.gz
cp $SUBJECTS_DIR/${visit}${subject}/vFREY.feat/reg/highres2standard_warp.nii.gz $SUBJECTS_DIR/${visit}${subject}/ASL_${cond}/highres2standard_warp.nii.gz
cp $SUBJECTS_DIR/${visit}${subject}/vFREY.feat/reg/standard2highres.mat $SUBJECTS_DIR/${visit}${subject}/ASL_${cond}/standard2highres.mat
cp $SUBJECTS_DIR/${visit}${subject}/vFREY.feat/reg/highres2standard.mat $SUBJECTS_DIR/${visit}${subject}/ASL_${cond}/highres2standard.mat


cd $SUBJECTS_DIR/${visit}${subject}/ASL_${cond}

# create transform between highres and ASL space

flirt -in calibration_body.nii.gz -ref ../strucfiles.anat/T1_biascorr.nii.gz -out rawASL${cond}2STRUC.nii.gz -omat rawASL${cond}2STRUC.mat -bins 256 -cost corratio -searchrx -30 30 -searchry -30 30 -searchrz -30 30 -dof 6  -paddingsize 0 -interp trilinear

convert_xfm -omat STRUC2rawASL${cond}.mat -inverse rawASL${cond}2STRUC.mat

### FLIRT brain-extracted T1 to ASL space 
flirt -in ../strucfiles.anat/T1_biascorr_brain.nii.gz -applyxfm -init STRUC2rawASL${cond}.mat -out rawASL${cond}BET.nii.gz -paddingsize 0 -interp trilinear -ref calibration_body.nii.gzs

### remove non brain tissue from magnitude file by masking with BETed structural images because it works better than BET
### flirt magnitude image to T1 nd then reverse it
flirt -in mag_ASL_${cond}.nii.gz -ref ../strucfiles.anat/T1_biascorr.nii.gz -out magASL${cond}2STRUC.nii.gz -omat magASL${cond}2STRUC.mat -bins 256 -cost corratio -searchrx -30 30 -searchry -30 30 -searchrz -30 30 -dof 6  -paddingsize 0 -interp trilinear

convert_xfm -omat STRUC2magASL${cond}.mat -inverse magASL${cond}2STRUC.mat

flirt -in ../strucfiles.anat/T1_biascorr_brain.nii.gz -applyxfm -init STRUC2magASL${cond}.mat -out magASL${cond}BET.nii.gz -paddingsize 0 -interp trilinear -ref mag_ASL_${cond}.nii.gz

fslmaths magASL${cond}BET.nii.gz -bin magASL${cond}BETBIN.nii.gz
 
fslmaths mag_ASL_${cond}.nii.gz -mas magASL${cond}BETBIN.nii.gz mag_ASL_${cond}_brain.nii.gz
 
### prepare FIELDMAP
fsl_prepare_fieldmap SIEMENS pha_ASL_${cond} mag_ASL_${cond}_brain.nii.gz fmap_ASL_${cond}.nii.gz 2.46

###  asl_mask1 is the input for GMfromPERFUSION1.m
### added 1 to names to tell the first, preliminary realignement to uncorrected ASL space from transforms to the corrected ASL space
fslmaths rawASL${cond}BET.nii.gz -bin asl_mask1.nii.gz

### create good anatomical contrast image in ASL space before fmap and bias correction:INPUTS: asl_raw_data and asl_mask1; OUTPUT: data1_std.nii.gz
matlab < ~karolina/Documents/MATLAB/GMfromPERFUSION1.m

### give the epi_reg the segmented white matter and run epi_reg
fslmaths ../strucfiles.anat/T1_fast_pve_2.nii.gz -thr 0.5 -bin struc_brain_wmseg.nii.gz

imcp struc_brain_wmseg.nii.gz asl2struc_fast_wmseg.nii.gz 

epi_reg --epi=data_std1.nii.gz --t1=../strucfiles.anat/T1_biascorr.nii.gz --t1brain=../strucfiles.anat/T1_biascorr_brain.nii.gz --out=asl2struc1 --fmap=fmap_ASL_${cond}.nii.gz --fmapmag=mag_ASL_${cond}.nii.gz --fmapmagbrain=mag_ASL_${cond}_brain.nii.gz --echospacing=0.00054 --pedir=-y --wmseg=asl2struc_fast_wmseg.nii.gz


### create a mask in ASL space after fieldmap correction
convert_xfm -omat struc1_2asl.mat -inverse asl2struc1.mat

flirt -in ../strucfiles.anat/T1_biascorr_brain.nii.gz -applyxfm -init struc1_2asl.mat -out asl_mask2NONbin.nii.gz -paddingsize 0 -interp trilinear -ref data_std1.nii.gz

fslmaths asl_mask2NONbin.nii.gz -bin asl_mask2.nii.gz

### bias-field correct the calibration files and the asl_raw_data and register calibration files to asl_raw_data

applywarp --postmat=asl2struc1_inv.mat --warp=asl2struc1_warp --in=asl_raw_data --ref=data_std1 --out=asl_raw_data_corrected_native

applywarp --postmat=asl2struc1_inv.mat --warp=asl2struc1_warp --in=calibration_head.nii.gz --ref=data_std1 --out=calibration_head_corrected_native

applywarp --postmat=asl2struc1_inv.mat --warp=asl2struc1_warp --in=calibration_body.nii.gz --ref=data_std1 --out=calibration_body_corrected_native

### MOTION CORRECTION of ASL data - has to be before epi_reg because asl_raw_data_corrected_mcf is an input for matlab
mcflirt -in asl_raw_data_corrected_native -out asl_raw_data_corrected_mcf -report

### create good anatomical contrast image in ASL space after fmap and bias correction (corrected ASL space); INPUTS: asl_raw_data_corrected_mcf and asl_mask2; OUTPUT: data_std.nii.gz: data_std.nii.gz
matlab < ~karolina/Documents/MATLAB/GMfromPERFUSION2.m

### The final transform between corrected ASL data and structural space
epi_reg --epi=data_std.nii.gz --t1=../strucfiles.anat/T1_biascorr.nii.gz --t1brain=../strucfiles.anat/T1_biascorr_brain.nii.gz --out=asl2struc --wmseg=asl2struc_fast_wmseg.nii.gz

convert_xfm -omat struc2asl.mat -inverse asl2struc.mat 

flirt -in ../strucfiles.anat/T1_biascorr_brain -applyxfm -init struc2asl.mat -out asl_maskNONbin.nii.gz -paddingsize 0 -interp trilinear -ref data_std.nii.gz

fslmaths asl_maskNONbin.nii.gz -bin asl_mask.nii.gz


### SUBTRACTION and AVERAGING of ASL; bias field correction removed - it is replaced by sensitivity correction later on 
asl_file --data=asl_raw_data_corrected_mcf --ntis=6 --ibf=rpt --iaf=tc --diff --out=asl_diff --mean=asl_diff_mean --mask=asl_mask.nii.gz

# threshold to remove negative voxels - mainly noise
fslmaths asl_diff_mean -thr 0 asl_diff_mean_thr


### Run asl_calib
### TR changed from 6 to 3.4 because of a sequence bug
### MAC: With the bias field correction above I think that the sensitivity correction is no longer called for (remove --cref and --osen):
### sensitivity correction added
asl_calib \
	-c calibration_head_corrected_native.nii.gz \
	-s struc_brain \
	-t asl2struc.mat \
	-m csfmask_asl \
	--mode longtr \
	--tissref csf \
	--te 0.013 \
	-o output \
	--t1r 4.3 \
	--t2r 0.75 \
	--tr 3.4 \
	--cgain 1 \
	--str2std highres2standard.mat \
	--warp highres2standard_warp.nii.gz \
	--cref calibration_body_corrected_native.nii.gz \
	--osen output/sens \
	--Mo output/M0.txt

SENS=output/sens
M0=`cat output/M0.txt`
ALPHA=0.88

### Correct perfusion data for coil sensitivity, M0 and inversion efficiency and convert into physiological units (ml/100/min) 
fslmaths asl_diff_mean_thr -div $SENS -mas $SENS -div $M0 -div $ALPHA output/perf_calib

### Transform the white and grey matter masks to the corrected ASL space  April 2015: changed spline to trilinear 
applywarp --ref=data_std.nii.gz --in=../strucfiles.anat/T1_fast_pve_1.nii.gz --out=pvgm_inasl_unthr --premat=struc2asl.mat --super --interp=spline --superlevel=4

applywarp --ref=data_std.nii.gz --in=../strucfiles.anat/T1_fast_pve_2.nii.gz --out=pvwm_inasl_unthr --premat=struc2asl.mat --super --interp=spline --superlevel=4

### threshold (upper and lower) the PVE to avoid artefacts of spline interpolation and also ignore very low PVE that could cause numerical issues. Threshold changed from 0.1 to 0.8 after MAC email Sept 2014
    fslmaths pvgm_inasl_unthr -thr 0.25 -min 1 pvgm_inasl
    fslmaths pvwm_inasl_unthr -thr 0.25 -min 1 pvwm_inasl

                
oxford_asl -i asl_diff_mean_thr \
		--tis 1.65,1.9,2.15,2.4,2.65,2.9 \
		-o output_vox \
		-s struc_brain.nii.gz \
        --asl2struc asl2struc.mat \
		--regfrom data_std \
		-t highres2standard.mat \
		--structout \
		--vars \
		--spatial \
		--casl \
		--bolus 1.4 \
		--slicedt 0.0452 \
		--fixbolus \
		--report \
		-m asl_mask \
		--bat 1.3 \
		--alpha 0.88 \
		--te 0.013 \
		--cgain 1 \
        --cmethod voxel \
        -c calibration_head_corrected_native.nii.gz \
		--pvgm pvgm_inasl \
		--pvwm pvwm_inasl \


### Convert perfusion maps into physiological units (ml/100/g/min)  6000 is to convert ml/g/s to ml/100g/min


fslmaths output/native_space/perfusion -mul 6000 output/native_space/perfusion
fslmaths output/struct_space/perfusion -mul 6000 output/struct_space/perfusion

fslmaths output/native_space/perfusion_var -mul 36000000 output/native_space/perfusion_var
fslmaths output/struct_space/perfusion_var -mul 36000000 output/struct_space/perfusion_var





done
done
done



