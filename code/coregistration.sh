
 dir=/media/hello/7T_Data

 for sub in `cat ${dir}/subject.list`
 do
	
 	cd ${dir}/${sub}
 	echo "${sub}"

 N4BiasFieldCorrection -d 3 -i b0_brain.nii.gz -o b0_brain_n4.nii.gz
 DenoiseImage -d 3 -i b0_brain_n4.nii.gz -o b0_brain_n4_de.nii.gz
	
 cp -f b0_brain_n4_de.nii.gz ${dir}/tmp/${sub}_b0.nii.gz
    
 done

   cd ${dir}/template
   antsMultivariateTemplateConstruction2.sh -d 3 -g 0.1 -o B0_temp -r 1 *.nii.gz

dir=/media/hello/data1

for sub in `cat ${dir}/subject.list`
do
	cd ${dir}/${sub}
	# 3dcopy ADC.hdr ADC.nii

	 antsApplyTransforms -d 3 -e 0 -i ADC.nii -r ${dir}/tmp/B0_temptemplate0.nii.gz -o ADC_warped.nii -t B0_temp*1Warp.nii.gz -t B0_temp*0GenericAffine.mat
done 

 dir=/media/hello//ROI

 for roi in `cat ${dir}/roi.list`
 do
   cd ${dir}
   antsApplyTransforms -d 3 -e 0 -i ${roi} -r B0_temptemplate0.nii.gz -o ${roi}_warped.nii -t TMBTA2tmp1Warp.nii.gz -t TMBTA2tmp0GenericAffine.mat
   fslmaths ${roi}_warped.nii -bin ${roi}_warped_bin.nii
 done

# ls >> .list

    dir=/media/hello/7T_Data/Ttest
	mkdir aqp_con
	3dttest++ -setA aqp \ -setB con \ -prefix ${dir}/aqp_con/aqp_con.nii.gz \ -mask ${dir}/B0_temptemplate0_mask.nii.gz