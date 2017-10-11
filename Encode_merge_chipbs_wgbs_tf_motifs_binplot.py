
import numpy as np
import pandas as pd
import pybedtools
import re
import os
import shutil
import time
import sys
import errno
from glob import glob
from pyfasta import Fasta
from os.path import basename
from os.path import join
from os.path import splitext

start_time = time.time()


### Can be chromHMM file or IDEAS segmentation file:
ideas_hepg2_file = os.path.expanduser("~/for_chris/batch_I/hepg2_ideas_36_dense.bed")

# chromHMM_hepg2_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/chrom_impute_chromHMM_data/corrected/E118_25_imputed12marks_dense.bed")
# refgene_file = os.path.expanduser("/Users/suryachhetri/Dropbox/local_miscellaneous_data/hg_19_misc/refGene_hg19")
# fasta_file = os.path.expanduser("~/for_chris/hg19-male.fa")
# fasta_idx = Fasta(fasta_file)

cpg_file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/merged_chipbs_rep1_rep2_w_wgbs" 
cpg_file_list = glob(join(cpg_file_dir, "*.bed"))
cpg_file_dict = {basename(each).split("_")[0]:each for each in cpg_file_list}

output_dir = os.path.expanduser("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/chipbs_tf_intersect")
final_output_file_1 = "_chip_wgbs_motifs_intersect_20bp_groupedby_bins.bed"
final_output_file_2 = "_chip_wgbs_motifs_intersect_20bp_each_cpgs.bed"

sub_output_dir_1 = join(output_dir, "binned_motifs_coord_info")
if not os.path.exists(sub_output_dir_1):
    os.makedirs(sub_output_dir_1)

sub_output_dir_2 = join(output_dir, "sorted_common_cpgs_info")
if not os.path.exists(sub_output_dir_2):
    os.makedirs(sub_output_dir_2)

sub_output_dir_3 = join(output_dir, "cpgs_and_tfs_binned_motif_intersect")
if not os.path.exists(sub_output_dir_3):
    os.makedirs(sub_output_dir_3)

### Needed for multiple threads running with a race condition to create the dir:
try:
    os.makedirs(output_dir)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass


def final_motifs_model(motif_file):
    motif_df = pd.read_csv(motif_file, sep="\t", header=None)
    ### Assuming that fimo motif and peaks 4th column has strand information:
    motif_select_df = motif_df.iloc[:, [0,1,2,3]]
    motif_select_df.columns = ["chrom", "start", "end", "strand"]
    print "\nCurrent dimension of the motif model:\n", motif_select_df.shape
    motif_select_df = motif_select_df.drop_duplicates()
    print "Dropping duplicates if any, current dimension of the motif model:\n", motif_select_df.shape
    return(motif_select_df)

motifs_coord_df = final_motifs_model("JunD_motifs_from_cpg_assoc_peaks.bed")


def generate_motifs_binned_coords(motifs_coordinates_info, tf_name, upstream_range, downstream_range, bin_size):
    motifs_df =  motifs_coordinates_info.sort_values(["chrom","start","end"])
    upstream = upstream_range
    downstream = downstream_range
    bin_size = bin_size
    nrows =  motifs_df.shape[0]

    bins = range(-upstream, (downstream), bin_size)
    bin_len = len(bins)
    motifs_concat_df = pd.concat([motifs_df]*bin_len, ignore_index="TRUE")
    motifs_sorted_df = motifs_concat_df.sort_values(["chrom","start","end"])

    ### Copy the bin list that is deep copy:
    bin_start_list = bins[:]
    bin_end_list = []
    for each in bin_start_list:
        bin_end_list.append(each+bin_size)

    bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
    bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
    bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

    ### Combine the motifs df and bin df by cbind or column wise:
    temp_motifs_df = pd.concat([motifs_sorted_df.reset_index(), bin_concat_df], axis = 1)
    temp_motifs_df["motifs_midpoint"] = (temp_motifs_df["start"] + temp_motifs_df["end"])/2
    temp_motifs_df["motifs_midpoint"] = temp_motifs_df["motifs_midpoint"].round().astype(int)
    final_motifs_df = temp_motifs_df.loc[:,["chrom", "motifs_midpoint", "bin_start", "bin_end", "start", "end"]]

    """ 
    chrom_start = tss_midpt + (bin_start); 
    chrom_end = tss_midpt + (bin_end) """
    final_motifs_df["chrom_start"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_start"]
    final_motifs_df["chrom_end"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_end"]

    select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'motifs_midpoint']
    final_motifs_df = final_motifs_df.loc[:,select_cols]

    ### Handle the case where your coordinates(mid point) are smaller than the upstream range (-2000) indicated; else start coordinate < 0(i.e -ve)
    final_motifs_df = final_motifs_df.loc[final_motifs_df["chrom_start"] > 0, :] 
    final_motifs_df.to_csv(join(sub_output_dir_1, tf_name + "_motifs_coordinate_info.bed"), sep="\t", index=False, header=False)

    return(final_motifs_df)

motifs_midpoint_coord_df= generate_motifs_binned_coords(motifs_coord_df, tf_name, 20, 20, 1)


def load_motifs_coord_pybedtool_object(file_name_with_full_path): 
    each_file = file_name_with_full_path
    os.environ["output_dir"] = sub_output_dir_1
    os.environ["each_file"] = each_file
    sorted_motif_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_motif_file"] = sorted_motif_file
    CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3,$4,$5,"."}' $each_file | sort -k1,1 -k2,2n  > $output_dir/$sorted_motif_file'''
    os.system(CMD)   
  
    print "Generating bedtools object..."
    motif_bed = pybedtools.BedTool(join(sub_output_dir_1,sorted_motif_file))
    
    return(motif_bed)

motifs_bed_file = load_motifs_coord_pybedtool_object(join(sub_output_dir_1, tf_name + "_motifs_coordinate_info.bed"))
    

def load_cpg_pybedtool_object(file_name_with_full_path):
    print " \nProcessing Cpg bed file\n "
    cpg_df = pd.read_csv(file_name_with_full_path, sep="\t")
    cpg_df_select = cpg_df[~(cpg_df["pval"] < 0.05)]
    cpg_df_select = cpg_df_select.iloc[:,0:10] # Removing the pval columns(10th col) from the file
    
    cpg_df_select["chip_meth"] = cpg_df_select["chipbs_rep1_meth"] + cpg_df_select["chipbs_rep2_meth"]
    cpg_df_select["chip_unmeth"] = cpg_df_select["chipbs_rep1_unmeth"] + cpg_df_select["chipbs_rep2_unmeth"]
    select_cols = ["chrom", "start", "end", "wgbs_meth", "wgbs_unmeth", "strand", "chip_meth", "chip_unmeth"]
    cpg_df_select = cpg_df_select.loc[:, select_cols]

    cpg_df_select.to_csv(join(sub_output_dir_2, tf_name + "_assoc_meth.bed"), sep="\t", header=True, index=False)
    ### Assuming the file's been sorted, important for faster intersect
    cpg_bed_file = pybedtools.BedTool.from_dataframe(cpg_df_select)

    return(cpg_bed_file)

cpg_bed_file = load_cpg_pybedtool_object(cpg_file_dict.get(tf_name))


def generate_motifs_binned_perc_meth(tf_name, motifs_final_bedfile, meth_file_list, **kwargs):
    # file_name =  kwargs["files_basename"]
    print "kwargs: ", kwargs
    meth_file_list = [cpg_file_dict.get(tf_name)]
    motifs_bedfile = motifs_bed_file

    master_dict = {}

    for idx, each_file in enumerate(meth_file_list):
        file_basename = splitext(basename(each_file))[0]    
        cpg_bed_file = load_cpg_pybedtool_object(each_file)

        """Pybedtool intersection of Cpg_bed_file and peak_bed_file"""
        print " Processing the intersection for", file_basename
        pybed_outfile = join(sub_output_dir_3, (tf_name + "_binned_motif_cpg_intersect.bed" ))

        #if not os.path.exists(pybed_outfile): 
        cpg_bed_intersect = cpg_bed_file.intersect(motifs_bedfile, wa = True, wb = True, output=pybed_outfile)      
        cpg_bed_intersect.head() 
        #else:
        #    print "Pybedtools object already present"

        if os.stat(pybed_outfile).st_size == 0: # handles the empty file or interrupted and killed file on prior jobs
            cpg_bed_intersect = cpg_bed_file.intersect(motifs_bedfile, wa = True, wb = True, output=pybed_outfile)      
            cpg_bed_intersect.head()
        else:
            print "Prior incomplete or undone intersection completed"

        ### Working with the dataframes; reading the output file of pybedtool intersect:
        df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
        df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12]]
        df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "wgbs_meth", "wgbs_unmeth", "chip_meth", "chip_unmeth", "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end"]    
        df_ordered[["wgbs_meth", "wgbs_unmeth", "chip_meth", "chip_unmeth"]] = df_ordered[["wgbs_meth", "wgbs_unmeth", "chip_meth", "chip_unmeth"]].astype(int)
        df_ordered["count"] = 1
        df_grouped =  df_ordered.groupby(["bin_start", "bin_end"]).apply(lambda x : (x["chip_meth"].sum()/float(x["chip_meth"].sum() + x["chip_unmeth"].sum()),\
                                                                                                                 x["wgbs_meth"].sum()/float(x["wgbs_meth"].sum() + x["wgbs_unmeth"].sum()),
                                                                                                                 x["count"].sum()
                                                                                                                 ))   

        df_grouped_reset = df_grouped.reset_index(name="percent_info")
        df_grouped_reset[["chip_meth_perc","wgbs_meth_perc", "cpg_count_in_bin"]] = df_grouped_reset["percent_info"].apply(pd.Series)

        print "Dimension of currently intersected peak file is", df_grouped.reset_index().shape

        if file_basename not in master_dict:
          master_dict[file_basename] = df_grouped_reset

        print "\nIntersection of %s completed!!!...\n\n" %(file_basename)

    ### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
    cpg_intersect_combined_df = pd.concat(master_dict).reset_index(drop=True)
    cpg_intersect_combined_df["annotation"] = tf_name + "_sites_meth_profile"
    cpg_intersect_combined_df.to_csv(join(output_dir, tf_name + final_output_file_1), sep ="\t", header = True, index = True)

    return(cpg_intersect_combined_df)


def generate_motifs_binned_perc_meth_for_each_cpgs(tf_name, motifs_final_bedfile, meth_file_list, **kwargs):
    print "kwargs: ", kwargs
    meth_file_list = [cpg_file_dict.get(tf_name)]
    motifs_bedfile = motifs_bed_file

    master_dict = {}

    for idx, each_file in enumerate(meth_file_list):
        file_basename = splitext(basename(each_file))[0]    
        cpg_bed_file = load_cpg_pybedtool_object(each_file)

        """Pybedtool intersection of Cpg_bed_file and peak_bed_file"""
        print " Processing the intersection for", file_basename
        pybed_outfile = join(sub_output_dir_3, (tf_name + "_binned_motif_cpg_intersect.bed" ))

        #if not os.path.exists(pybed_outfile): 
        cpg_bed_intersect = cpg_bed_file.intersect(motifs_bedfile, wa = True, wb = True, output=pybed_outfile)      
        cpg_bed_intersect.head() 
        #else:
        #    print "Pybedtools object already present"

        if os.stat(pybed_outfile).st_size == 0: # handles the empty file or interrupted and killed file on prior jobs
            cpg_bed_intersect = cpg_bed_file.intersect(motifs_bedfile, wa = True, wb = True, output=pybed_outfile)      
            cpg_bed_intersect.head()
        else:
            print "Prior incomplete or undone intersection completed"

        ### Working with the dataframes; reading the output file of pybedtool intersect:
        df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
        df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12]]
        df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "wgbs_meth", "wgbs_unmeth", "chip_meth", "chip_unmeth", "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end"]    
        df_ordered[["wgbs_meth", "wgbs_unmeth", "chip_meth", "chip_unmeth"]] = df_ordered[["wgbs_meth", "wgbs_unmeth", "chip_meth", "chip_unmeth"]].astype(int)
        df_ordered["cpg_count"] = 1

        df_ordered["chip_meth_perc"] = df_ordered["chip_meth"]/(df_ordered["chip_meth"] + df_ordered["chip_unmeth"])
        df_ordered["wgbs_meth_perc"] = df_ordered["wgbs_meth"]/(df_ordered["wgbs_meth"] + df_ordered["wgbs_unmeth"])
        df_ordered = df_ordered.loc[:, ["cpg_chrom", "cpg_start", "cpg_end", "bin_start", "bin_end", "chip_meth_perc", "wgbs_meth_perc", "cpg_count"]] 

        if file_basename not in master_dict:
          master_dict[file_basename] = df_ordered

        print "\nIntersection of %s completed!!!...\n\n" %(file_basename)

    ### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
    cpg_intersect_combined_df = pd.concat(master_dict).reset_index(drop=True)
    cpg_intersect_combined_df["annotation"] = tf_name + "_sites_meth_profile"
    cpg_intersect_combined_df.to_csv(join(output_dir, tf_name + final_output_file_2), sep ="\t", header = True, index = True)

    return(cpg_intersect_combined_df)


#######################################################


""" Motif methylation analysis of all TFs, since above codes was designed for single TF """
motif_file_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/filtered_motifs_from_cpg_assoc_peaks/*cpg_assoc_peaks.bed"
motif_file_list = glob(motif_file_path)
motif_file_dict = {basename(each).split("_")[0]:each for each in motif_file_list}
tf_name_list = ["JunD", "Pol2", "TAF1", "YY1"]

### This is useful to do run all the files parallell in the cluster at a time:
### Usage: bsub -We -n 1 -o "./fimo_motif_bin_methplot.out $RUNPATH/Encode_motifs_fimo_bin_methplot_final_total.py fimo file
# motif_file_list = [sys.argv[1]]

master_tf_dict_1 = {}
master_tf_dict_2 = {}
for tf_name in tf_name_list:
    print "\nCurrently processing %s tf\n" %(tf_name)

    # motifs_coord_df = final_motifs_model("JunD_motifs_from_cpg_assoc_peaks.bed")
    motifs_coord_df = final_motifs_model(motif_file_dict.get(tf_name))
    motifs_midpoint_coord_df= generate_motifs_binned_coords(motifs_coord_df, tf_name, 20, 20, 1)       
    
    motifs_bed_file = motifs_bed_file = load_motifs_coord_pybedtool_object(join(sub_output_dir_1, tf_name + "_motifs_coordinate_info.bed"))
    list_of_cpg_files = [cpg_file_dict.get(tf_name)]

    ### 2 types of data generated: a) avg. methylation on each bins/groupby bins 
    ### b) avg. methylation on each cpgs based on bins for detailed boxplot dist:
    plot_data_set_goupedby_bins = generate_motifs_binned_perc_meth(tf_name, motifs_bed_file, list_of_cpg_files)    
    plot_data_set_each_cpgs = generate_motifs_binned_perc_meth_for_each_cpgs(tf_name, motifs_bed_file, list_of_cpg_files)

    if tf_name not in master_tf_dict_1:
            master_tf_dict_1[tf_name] = plot_data_set_goupedby_bins

    if tf_name not in master_tf_dict_2:
            master_tf_dict_2[tf_name] = plot_data_set_each_cpgs

master_tf_df_groupedby_bins = pd.concat(master_tf_dict_1).reset_index(drop=True)
master_tf_df_each_cpg_based = pd.concat(master_tf_dict_2).reset_index(drop=True)

""" Final containing information of all the TFs motif methylation profile """
master_tf_df_groupedby_bins.to_csv(join(output_dir, "all_tfs_" + final_output_file_1), sep ="\t", header = True, index = True)
master_tf_df_each_cpg_based.to_csv(join(output_dir, "all_tfs_" + final_output_file_2), sep ="\t", header = True, index = True)

### Comment above line code and uncomment the below, if you running parallel in the cluster since master_tf_df dataframe is output of the motifs fimo bin methplot script result:
#master_tf_df.to_csv(join(output_dir, (tf_name + "_" + final_output_file)), sep ="\t", header = True, index = True)

print "Task completed!!!"
print "\n\nTime taken to complete the analayis", time.time()-start_time







# plot_output_dir = join(output_dir,"composite_methplots") # os.path.expanduser("~/for_chris/batch_I/motifs_methplot/composite_methplots")
# if not os.path.exists(plot_output_dir):
#     os.makedirs(plot_output_dir)

# wgbs_motif_intersect_file_1 = join(output_dir, tf_name + final_output_file_1)
# wgbs_motif_intersect_file_2 = join(output_dir, tf_name + final_output_file_2)
# os.environ["wgbs_motif_file_1"] =  wgbs_motif_intersect_file_1
# os.environ["wgbs_motif_file_2"] =  wgbs_motif_intersect_file_2
# os.environ["tf_name"] = tf_name
# os.environ["plot_output_dir"] = plot_output_dir

# os.system('Rscript ./Encode_motifs_bin_methplot_final.R $wgbs_motif_file_1 $tf_name $plot_output_dir')
# os.system('Rscript ./Encode_motifs_bin_methplot_final.R $wgbs_motif_file_2 $tf_name $plot_output_dir')
# print "\nRunning of Rscript for plotting completed!!!....\n"
# print "\nCheck your plots in %s dir\n" %(plot_output_dir)