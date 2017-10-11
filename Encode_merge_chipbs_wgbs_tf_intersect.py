import os
import glob
import re
import pandas as pd
import pybedtools
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as smm
from os.path import splitext
from os.path import basename
from os.path import join


### Pre-processing of data in bash:
# SL_list="SL60584 SL94837 SL68320 SL88936 SL81122 SL88934 SL83433 SL88937 SL68322"
# for each in $SL_list; do cp /gpfs/gpfs1/home/schhetri/final_chipBS_data/$each/unsortedButMerged_ForBismark_file/methylation_extraction/*.cov . ; done
# for each in $(ls *.cov); do file_name=$(basename $each .cov); sort -k1,1 -k2,2n $each > ${file_name}_sorted.cov; done
# cut -f2- seq-specific factors | sort -k1,1 -k2,2n > seq-specific-factors.bed (# Remove motif 1 notation from seq-sp factors)

tf_name_list = ["JunD", "YY1", "ZBTB33", "TAF1", "Pol2"]
tf_name_list = ["JunD","YY1", "TAF1", "Pol2"] # ZBTB has single replicate
# tf_name_list = ["JunD"] # ZBTB has single replicate

def meth_percent_calc(x):
    x_meth = x["cpg_meth"].sum()
    x_unmeth = x["cpg_meth"].sum() + x["cpg_unmeth"].sum()
    x_perc = x_meth/float(x_unmeth)
    return(x_perc)

### Two sided fisher test:
def fishers_test(df_rows):
    x1 = df_rows["chipbs_rep1_meth"]
    x2 = df_rows["chipbs_rep1_unmeth"]
    y1 = df_rows["chipbs_rep2_meth"]
    y2 = df_rows["chipbs_rep2_unmeth"]
    pvalue = stats.fisher_exact([[x1,x2],[y1,y2]])[1]

    """ #total_pval_test = df_rows["sig_hits"].shape[0]
        #bonferroni_correction = min(pvalue*total_pval_test, 1)
        #return(pvalue, df_rows["tf_name"]) """

    return pd.Series({"pvalue": pvalue})


### Input and output file and dir paths:
main_dir = os.path.expanduser("/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis")
if not os.path.exists(main_dir):
    os.makedirs(main_dir)

sorted_peak_file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/peak_files"
if not os.path.exists(sorted_peak_file_dir):
    os.makedirs(sorted_peak_file_dir)

peak_file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/raw_peak_files"
if not os.path.exists(peak_file_dir):
    os.makedirs(peak_file_dir)


peak_bed_list = glob.glob(join(peak_file_dir, "*.bed"))
peaks_file_dict = {basename(each).split("_")[0]:each for each in peak_bed_list}

cpg_assoc_peaks = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/cpg_associated_peaks"
if not os.path.exists(cpg_assoc_peaks):
    os.makedirs(cpg_assoc_peaks)

### Already created using sorted peak beds for non-sequence specfic factors and motif 1 for seq-specific factors:
### cut -f2- seq-specific factors | sort -k1,1 -k2,2n > seq-specific-factors.bed (# Remove motif 1 notation from seq-sp factors)
motifs_info = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/motifs_info"
motifs_info_file_list = glob.glob(join(motifs_info, "*bed")) 
motifs_info_file_dict = {basename(each).split("_")[0]:each for each in motifs_info_file_list}

filtered_motifs_from_cpg_assoc_peaks ="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/filtered_motifs_from_cpg_assoc_peaks" 
if not os.path.exists(filtered_motifs_from_cpg_assoc_peaks):
    os.makedirs(filtered_motifs_from_cpg_assoc_peaks)

chipbs_tf_dir = join(main_dir, "chipbs_tf_intersect")
if not os.path.exists(chipbs_tf_dir):
    os.makedirs(chipbs_tf_dir)

wgbs_tf_dir = join(main_dir, "merged_chipbs_rep_w_wgbs_tf_intersect")
if not os.path.exists(wgbs_tf_dir):
    os.makedirs(wgbs_tf_dir)

wgbs_file_dir =  join(main_dir, "wgbs_file")
if not os.path.exists(wgbs_file_dir):
    os.makedirs(wgbs_dir)

cpg_bed_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL83768_SL83771_300x_merged_CpG_context_deduplicated.bismark.cov"

""" Creating wgbs cpg pybedtool object """
os.environ["cpg_bed_file"] = cpg_bed_file
cpg_bed_edited = "wgbs_" + splitext(basename(cpg_bed_file))[0] + "_edited" + splitext(cpg_bed_file)[1]
os.environ["cpg_bed_edited"] = cpg_bed_edited
os.environ["output_dir"] = wgbs_file_dir

if not os.path.exists(join(wgbs_file_dir,cpg_bed_edited)):
    os.system("ln -fs $cpg_bed_file $output_dir/$cpg_bed_edited")

""" Merge chipbs-replicates and wgbs common CpGs and Filter statistically significantly different CpGs from the replicates """
merged_chipbs_wgbs_file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/merged_chipbs_rep1_rep2_w_wgbs"

if not os.path.exists(merged_chipbs_wgbs_file_dir):
    os.makedirs(merged_chipbs_wgbs_file_dir)

chipbs_rep1 = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/chipbs_rep1"
chipbs_rep1_file_list = glob.glob(join(chipbs_rep1, "*cov")) 
chipbs_rep1_file_dict = {basename(each).split("_")[0]:each for each in chipbs_rep1_file_list}

chipbs_rep2 = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_bs_wgbs_analysis/chipbs_rep2"
chipbs_rep2_file_list = glob.glob(join(chipbs_rep2, "*cov")) 
chipbs_rep2_file_dict = {basename(each).split("_")[0]:each for each in chipbs_rep2_file_list}

wgbs_cpg_bed_file = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/SL83768_SL83771_300x_merged_CpG_context_deduplicated.bismark.cov"

master_tf_dict = {}
for each_tf in tf_name_list:
    TF_name = each_tf
    print "Processing : %s" %(TF_name)
   
    wgbs_cpg_df = pd.read_csv(wgbs_cpg_bed_file, sep = "\t", header = None)
    wgbs_cpg_df = wgbs_cpg_df.iloc[:, [0, 1, 2, 3, 4]]
    wgbs_cpg_df.columns = ["chrom", "start", "end", "wgbs_meth", "wgbs_unmeth"]
    
    chipbs_rep1_df = pd.read_csv(chipbs_rep1_file_dict.get(TF_name), sep = "\t", header = None)
    chipbs_rep1_df = chipbs_rep1_df.iloc[:, [0, 1, 2, 3, 4]]
    chipbs_rep1_df.columns = ["chrom", "start", "end", "chipbs_rep1_meth", "chipbs_rep1_unmeth"]

    chipbs_rep2_df = pd.read_csv(chipbs_rep2_file_dict.get(TF_name), sep = "\t", header = None)
    chipbs_rep2_df = chipbs_rep2_df.iloc[:, [0, 1, 2, 3, 4]]
    chipbs_rep2_df.columns = ["chrom", "start", "end", "chipbs_rep2_meth", "chipbs_rep2_unmeth"]

    all_dfs = [wgbs_cpg_df, chipbs_rep1_df, chipbs_rep2_df]

    if not os.path.exists(join(merged_chipbs_wgbs_file_dir, TF_name + "_chipbs_wgbs_merged_cpgs.bed")):
        merged_cpg_df = reduce(lambda left, right: pd.merge(left, right, on=["chrom","start","end"]), all_dfs)
        pval_df = merged_cpg_df.apply(fishers_test, axis=1)
        merged_cpg_df["pval"] = pval_df["pvalue"]
        merged_cpg_df["strand"] = "."
        merged_cpg_df_select = merged_cpg_df.loc[:,[u'chrom', u'start', u'end', u'wgbs_meth', u'wgbs_unmeth', u'strand',\
                                    u'chipbs_rep1_meth', u'chipbs_rep1_unmeth', u'chipbs_rep2_meth', u'chipbs_rep2_unmeth', u'pval']]
        
        merged_cpg_df_select.to_csv(join(merged_chipbs_wgbs_file_dir, TF_name + "_chipbs_wgbs_merged_cpgs.bed"), sep ="\t", header = True, index = False)

        # pvalue_list = pval_df["pvalue"].tolist()
        # rej, pval_corr = smm.multipletests(pvalue_list, alpha=0.05, method="fdr_bh")[:2] # method="bonferroni" or "hs"; if needed
        # pval_df["pval_corr"] = list(pval_corr)
    
    ### Filter significantly different CpGs with alpha=0.025(2 tail test) with confidence level of 95% :
    merged_chipbs_wgbs_df = pd.read_csv(join(merged_chipbs_wgbs_file_dir, TF_name + "_chipbs_wgbs_merged_cpgs.bed"), sep="\t")
    merged_chipbs_wgbs_df_select = merged_chipbs_wgbs_df[~(merged_chipbs_wgbs_df["pval"] < 0.05)]
    merged_chipbs_wgbs_file = pybedtools.BedTool.from_dataframe(merged_chipbs_wgbs_df_select)
    merged_chipbs_wgbs_file.head()

    peaks_file_path = peaks_file_dict.get(TF_name)
    print "\nPeak file path for %s : %s\n" %(TF_name, peaks_file_path)

    """ \nProcessing and creating a sorted bed file"""  
    os.environ["each_file"] = peaks_file_path
    os.environ["output_dir"] = sorted_peak_file_dir
    sorted_peak_file = splitext(basename(peaks_file_path))[0] + "_sorted" + splitext(basename(peaks_file_path))[1]
    os.environ["sorted_peak_file"] = sorted_peak_file
    if not os.path.exists(join(sorted_peak_file_dir, sorted_peak_file)):
        os.system('sort -k1,1 -k2,2n $each_file > $output_dir/$sorted_peak_file')

    peaks_bed_file = pybedtools.BedTool(join(sorted_peak_file_dir, sorted_peak_file))
    print peaks_bed_file.head()

    """ \nPybedtool intersection of Cpg_bed_file and peak_bed_file...\n """
    print " Processing the intersection for", TF_name
    pybed_outfile = join(wgbs_tf_dir, (TF_name + "_sorted_cpg_intersect_wgbs.bed"))
    if not os.path.exists(pybed_outfile):
        wgbs_tf_intersect = merged_chipbs_wgbs_file.intersect(peaks_bed_file, wa = True, wb = True, output=pybed_outfile)
        print wgbs_tf_intersect.head()
    else:
        print "Pybedtools object already present"

    ### Working with the dataframes; reading the output file of pybedtool intersect:
    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
    df_ordered = df_file.iloc[:,[11,12,13,0,1,2,3,4,6,7,8,9]]
    df_ordered.columns = ["peak_chrom", "peak_start", "peak_end", "cpg_chrom", "cpg_start", "cpg_end", "wgbs_meth", "wgbs_unmeth", "chip1_meth", "chip1_unmeth", "chip2_meth", "chip2_unmeth"]
    
    df_ordered["chip_meth"] = df_ordered["chip1_meth"] + df_ordered["chip2_meth"]
    df_ordered["chip_unmeth"] = df_ordered["chip1_unmeth"] + df_ordered["chip2_unmeth"]
    df_ordered["chip_coverage"] = df_ordered["chip_meth"] + df_ordered["chip_unmeth"]
    df_ordered["wgbs_coverage"] = df_ordered["wgbs_meth"] + df_ordered["wgbs_unmeth"]
    df_ordered["count"] = 1

    ### Grouping each peaks to find weighted mean methylation and chipbs-wgbs ratio, basically peak wise calculation:
    df_grouped = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x : (x["chip_meth"].sum()/float(x["chip_meth"].sum() + x["chip_unmeth"].sum()),\
                                                                                                             x["wgbs_meth"].sum()/float(x["wgbs_meth"].sum() + x["wgbs_unmeth"].sum()),
                                                                                                             x["chip_coverage"].sum(),
                                                                                                             x["wgbs_coverage"].sum(),
                                                                                                             x["count"].sum()
                                                                                                             ))   
    df_grouped_reset = df_grouped.reset_index(name="percent_info")
    df_grouped_reset[["chip_meth_perc","wgbs_meth_perc", "total_chip_coverage", "total_wgbs_coverage", "cpg_count"]] = df_grouped_reset["percent_info"].apply(pd.Series)
    df_grouped_reset["chip_wgbs_ratio"] = df_grouped_reset["chip_meth_perc"]/df_grouped_reset["wgbs_meth_perc"]
    df_grouped_reset["cpg_count"] = df_grouped_reset["cpg_count"].astype(int)

    ### Normalization for the wgbs ratio (round each fraction to 2 decimal places and replace 0 by 0.001):
    df_grouped_reset["chip_meth_perc_norm"] = df_grouped_reset["chip_meth_perc"].round(2).replace({0:0.001})
    df_grouped_reset["wgbs_meth_perc_norm"] = df_grouped_reset["wgbs_meth_perc"].round(2).replace({0:0.001})
    df_grouped_reset["chip_wgbs_ratio_norm_log"] = np.log10(df_grouped_reset["chip_meth_perc_norm"]/df_grouped_reset["wgbs_meth_perc_norm"])
    df_grouped_reset.to_csv(join(chipbs_tf_dir, TF_name + "_peak_grouped_chipbs_wgbs_ratio.txt"), sep="\t", header=True, index=False)

    ### Calculate mean methylation for each cpgs overlapping the peaks, basically cpg wise calculation:
    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
    cpg_ordered = df_file.iloc[:,[0,1,2,3,4,6,7,8,9]]
    cpg_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "wgbs_meth", "wgbs_unmeth", "chip1_meth", "chip1_unmeth", "chip2_meth", "chip2_unmeth"]
    
    cpg_ordered["chip_meth"] = cpg_ordered["chip1_meth"] + cpg_ordered["chip2_meth"]
    cpg_ordered["chip_unmeth"] = cpg_ordered["chip1_unmeth"] + cpg_ordered["chip2_unmeth"]
    cpg_ordered["chip_coverage"] = cpg_ordered["chip_meth"] + cpg_ordered["chip_unmeth"]
    cpg_ordered["wgbs_coverage"] = cpg_ordered["wgbs_meth"] + cpg_ordered["wgbs_unmeth"]
    cpg_ordered.drop_duplicates(["cpg_chrom", "cpg_start", "cpg_end"], inplace=True)

    cpg_ordered["chip_meth_perc"] = cpg_ordered["chip_meth"]/(cpg_ordered["chip_meth"] + cpg_ordered["chip_unmeth"])
    cpg_ordered["wgbs_meth_perc"] = cpg_ordered["wgbs_meth"]/(cpg_ordered["wgbs_meth"] + cpg_ordered["wgbs_unmeth"])
    cpg_ordered["chip_wgbs_ratio"] = cpg_ordered["chip_meth_perc"]/cpg_ordered["wgbs_meth_perc"]

    ### Normalization for the wgbs ratio (round each fraction to 2 decimal places and replace 0 by 0.001):
    cpg_ordered["chip_meth_perc_norm"] = cpg_ordered["chip_meth_perc"].round(2).replace({0:0.001})
    cpg_ordered["wgbs_meth_perc_norm"] = cpg_ordered["wgbs_meth_perc"].round(2).replace({0:0.001})
    cpg_ordered["chip_wgbs_ratio_norm_log"] = np.log10(cpg_ordered["chip_meth_perc_norm"]/cpg_ordered["wgbs_meth_perc_norm"])
    cpg_ordered = cpg_ordered.loc[:, ["cpg_chrom", "cpg_start", "cpg_end", "chip_meth_perc", "chip_meth_perc_norm", "wgbs_meth_perc", 
                                        "wgbs_meth_perc_norm", "chip_wgbs_ratio", "chip_wgbs_ratio_norm_log", "chip_coverage","wgbs_coverage"]]
    cpg_ordered.to_csv(join(chipbs_tf_dir, TF_name + "_cpg_grouped_chipbs_wgbs_ratio.txt"), sep="\t", header=True, index=False)

    ### Grouping the cpg and save to the dir for further intersection with motifs for motif based analysis using chipbs-wgbs data:
    df_grouped_peaks = df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).size().reset_index(name="peak_count")
    df_grouped_peaks.to_csv(join(cpg_assoc_peaks, TF_name + "_cpg_associated_peaks.bed"), sep="\t", header=True, index=False)

    ### Read the cpg_associated_peaks file :
    cpg_peak_df = pd.read_csv(join(cpg_assoc_peaks, TF_name + "_cpg_associated_peaks.bed"), sep="\t")
    cpg_peak_df = cpg_peak_df.iloc[:,[0,1,2]]

    ### Filter the motifs or peaks overlapping with the cpg associated peaks only:
    cpg_peak_bedfile = pybedtools.BedTool.from_dataframe(cpg_peak_df)  
    motifs_info_bedfile = pybedtools.BedTool(motifs_info_file_dict.get(TF_name))
    pybed_outfile = join(filtered_motifs_from_cpg_assoc_peaks, (TF_name + "_motifs_from_cpg_assoc_peaks.bed"))
    if not os.path.exists(pybed_outfile):
        motif_cpg_peak_intersect = motifs_info_bedfile.intersect(cpg_peak_bedfile, u=True, output=pybed_outfile)
        print motif_cpg_peak_intersect.head()

print "Job completed successfully!!!"



# if "fimo" not in motifs_info_file_dict.get(TF_name):
# peak_df = pd.read_csv(motifs_info_file_dict.get(TF_name), sep = "\t", header = None)
#     peak_df = peaks_df.iloc[:,[0,1,2]]
#     peak_df.columns = ["chrom", "start", "end"]
#     peak_df["midpoint"] = (peak_df["start"] + peak["end"])/2
#     peak_df["midpoint"] = peak_df["midpoint"].round().astype(int)
# motifs_info_df = pd.read_csv(motifs_info_file_dict.get(TF_name), sep = "\t", header = None)

# df_ordered["chip1_meth_perc"] = df_ordered["chip1_meth"]/(df_ordered["chip1_meth"] + df_ordered["chip1_unmeth"])
# df_ordered["chip2_meth_perc"] = df_ordered["chip2_meth"]/(df_ordered["chip2_meth"] + df_ordered["chip2_unmeth"])
# df_ordered["wgbs_meth_perc"] = df_ordered["wgbs_meth"]/(df_ordered["wgbs_meth"] + df_ordered["wgbs_unmeth"])
# df_ordered["chip1_coverage"] = df_ordered["chip1_meth"] + df_ordered["chip1_unmeth"]
# df_ordered["chip2_coverage"] = df_ordered["chip2_meth"] + df_ordered["chip2_unmeth"]
# df_ordered["wgbs_coverage"] = df_ordered["wgbs_meth"] + df_ordered["wgbs_unmeth"]
# df_ordered["avg_chipbs_perc"] = (df_ordered["chip1_meth_perc"] + df_ordered["chip2_meth_perc"])/2
# df_grouped =  df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x : (x["wgbs_meth_perc"].mean(), x["avg_chipbs_perc"].mean()))
# df_grouped = df_grouped.reset_index(name="percent_info")
# df_grouped[["avg_wgbs_percent", "avg_chipbs_percent"]] = df_grouped["percent_info"].apply(pd.Series)
# df_grouped =  df_ordered.groupby(["peak_chrom", "peak_start", "peak_end"]).apply(lambda x : (x["wgbs_meth_perc"].mean(), x["avg_chipbs_perc"].mean()))
# df_grouped = df_grouped.reset_index(name="percent_info")
# df_grouped[["avg_wgbs_percent", "avg_chipbs_percent"]] = df_grouped["percent_info"].apply(pd.Series)

# df_grouped = df_ordered.groupby(["cpg_chrom", "cpg_start", "cpg_end"]).apply(lambda x : (x["chip_meth"].sum()/float(x["chip_meth"].sum() + x["chip_unmeth"].sum()),\
#                                                                                                          x["wgbs_meth"].sum()/float(x["wgbs_meth"].sum() + x["wgbs_unmeth"].sum()),
#                                                                                                          x["chip_coverage"].sum(),
#                                                                                                          x["wgbs_coverage"].sum(),
#                                                                                                          x["count"].sum()
#                                                                                                          ))   
# df_grouped_reset = df_grouped.reset_index(name="percent_info")
# df_grouped_reset[["chip_meth_perc","wgbs_meth_perc", "total_chip_coverage", "total_wgbs_coverage", "cpg_count"]] = df_grouped_reset["percent_info"].apply(pd.Series)
# df_grouped_reset["chip_wgbs_ratio"] = df_grouped_reset["chip_meth_perc"]/df_grouped_reset["wgbs_meth_perc"]
# df_grouped_reset["chip_wgbs_ratio"] = df_grouped_reset["chip_wgbs_ratio"].astype(str)
# df_grouped_reset["cpg_count"] = df_grouped_reset["cpg_count"].astype(int)


