##Bash file that will annotate and filter .vcf files of variants from sequenced genomes.
##There is another copy of this in the annovar folder on the desktop.
##Should be run from a linux shell within the annovar folder on the desktop.
##You will end up with .tsv files with two unnamed columns at the end, which are (left) read depth and (right) Variant Allele Frequency.
##Instructions to use this code:
##Put all patient vcf files in the folder "patient_vcf".
##In the annovar folder, ctrl+shift+right click in folder window, open linux shell, type command "bash annotate_automatically_LB.sh", hit enter.
##After a while the results should appear in the folder annovar/patient_vcf/patient_vcf_qc/patient_vcf_qc_anno
##This can take quite a while. If everything is working you should see files appearing and disappearing in patient_vcf_qc, and see
##things happening in the linux shell. You will know it is complete when the linux shell 
##goes back to the green text showing the address.
##Below is the code that will carry out the annotation and filtering.
#
##Step 1: Remove all entries where there is no base change (shown by the 'alternative' base in column 5 being ".")
#
cd patient_vcf
for FILE in *.vcf; 
do awk '$5!="."' $FILE > patient_vcf_qc/$FILE"_IDremoved.vcf"; 
done
#
##Step 2: quality control: remove poor reads from the file.
#
cd patient_vcf_qc
for FILE in *.vcf; 
do grep -v -e "Blacklist" -e "LowDP" -e "LowSupport" -e "LowVarSupport" -e "q0" $FILE > $FILE"_qc.vcf"; 
rm -f $FILE 
done
#
##Step 3: Annotate! Outputs to the same folder
#
cd ..
cd ..
for FILE in patient_vcf/patient_vcf_qc/*.vcf;  
do ./table_annovar.pl $FILE /mnt/c/Users/liam.brockley.UC/Desktop/annovar/humandb/ -buildver hg19 -out $FILE"_anno" -remove -protocol knownGene,gnomad211_genome,clinvar_20210501,cosmic70 -operation g,f,f,f -nastring . -vcfinput; 
done
cd patient_vcf/patient_vcf_qc
mv *.txt patient_vcf_qc_anno
#
##Extra: convert text to .tsv files with simple names
#
cd patient_vcf_qc_anno
for FILE in *.txt; 
do mv $FILE ${FILE%.vcf_IDremoved*}_anno.tsv; 
done
#
##The final column has a bunch of information separated by colons. This splits these all into new columns.
#
for FILE in *.tsv; 
do less $FILE | sed 's/ /_/g' | cut -f 46 | sed 's/:/\t/g' | paste $FILE - > $FILE"_runinfo.tsv"; 
rm $FILE 
done
#
##Now finally we remove unnecessary columns
##Columns 12-26 contain population allele frequencies in different demographics - we only need overall AF so I remove these
##Column 28 contains some ID number that we don't really need
##Columns 34-49 contain duplicate information
##Columns 52-57 contain unnecessary information from the run 
##This cutting just leaves you with the read depth and the VAF (for the patient), in the last 2 columns
#
for FILE in *.tsv; 
do cut --complement -f12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,52,53,54,55,56,57 $FILE > $FILE"_inforemoved.tsv";  
rm $FILE
done
#
##Now to rename everything again - a bit redundant really
#
for FILE in *.tsv; 
do mv $FILE ${FILE%_anno*}_.tsv; 
done
#
##There is still an annoying space that shifts columns over that must be removed in this step
#
for FILE in *.tsv; 
do sed 's/ /_/' $FILE > $FILE"_final.tsv"; 
rm $FILE
done
#
#Second for loop makes a second file for each that includes only variants with population frequency at or below 0.01, which will likely be somatic
#
for FILE in *.tsv; 
do awk '$11=="AF" || $11=="." || $11<=0.01' $FILE > $FILE"_lowPopAF.tsv"; 
done