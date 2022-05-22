#!/usr/bin/env python

import numpy as np
import pandas as pd
import pybedtools
import re, os, shutil
import time
import sys, errno
from glob import glob
from pyfasta import Fasta
from os.path import basename, join, splitext


start_time = time.time()

def final_fimo_motifs_model(fimo_file, tf_name, **kwargs):
    motif_list = []
    for key, values in kwargs.iteritems():
        motif_list.append(values)

    df = pd.read_csv(fimo_file, sep="\t")
    df.rename(columns={"sequence name" : "chrom", "#pattern name" : "motif_name"}, inplace=True)
    df = df.loc[df["motif_name"].isin(motif_list)]
    df["start"] = df["start"] #Extended 2 bp upstream
    df["end"] =  df["stop"] #Extended 2 bp downstream
    df["tf_name"] = tf_name
    df["motif_id"] = "MOTIF" + df["motif_name"].astype(str)
    select_cols = ["chrom", "start", "end", "motif_id", "tf_name", "p-value"]
    motif_select_df = df.loc[:,select_cols]
    motif_select_df = motif_select_df.loc[motif_select_df["p-value"] < 0.001]
    print("Current dimension of motif model : {}".format(motif_select_df.shape))
    # motif_select_df.duplicated(keep=False)
    motif_select_df = motif_select_df.drop_duplicates()
    print("Dropping duplicates if any, current dimension of motif model : {}\n".format(motif_select_df.shape))
    motif_sorted_df = motif_select_df.sort_values(["chrom", "start","end"]).reset_index(drop=True)
    return(motif_sorted_df)

# df =  final_fimo_motifs_model(fimo_file, "test",motif_1= 1, motif_2= 2)

def parse_fimo_motif_coords(fimo_file):
    df = pd.read_csv(fimo_file, sep="\t")
    master_dict = {}
    for index, row in df.iterrows():                                       
        motif_name = str(row["#pattern name"])
        chrom_name = row["sequence name"]
        start = row["start"]
        end = row["stop"]
        strand = row["strand"]
        fimo_seq = row["matched sequence"]

        if strand == "+":          
            start = start + 1
            end = end + 1
            seq_pyfasta = fasta_idx.sequence({"chr": chrom_name, "start" : start, "stop" : end, "strand" : strand}) # one based
                
        elif strand == "-":
            start = start + 1
            end = end + 1
            seq_pyfasta = fasta_idx.sequence({"chr": chrom_name, "start" : start, "stop" : end, "strand" : strand}) # one based

        ## print motif_name, chrom_name, start, end
        if motif_name not in master_dict:
            motif_name = motif_name
            master_dict[motif_name]={"chrom":[chrom_name],"start":[start], "end":[end], "strand": [strand], "fimo_seq": [fimo_seq], "pyfasta_seq":[seq_pyfasta]}
        else:
            master_dict[motif_name]["chrom"].append(chrom_name)
            master_dict[motif_name]["start"].append(start)
            master_dict[motif_name]["end"].append(end)
            master_dict[motif_name]["strand"].append(strand)
            master_dict[motif_name]["fimo_seq"].append(fimo_seq)
            master_dict[motif_name]["pyfasta_seq"].append(seq_pyfasta)
    return(master_dict)

# Master_motif_dict = parse_fimo_motif_coords(fimo_file)

def generate_motifs_binned_coords(motifs_coordinates_info, upstream_range, downstream_range, bin_size):
    motifs_df =  motifs_coordinates_info.sort_values(["chrom","start","end"])
    upstream = upstream_range
    downstream = downstream_range
    bin_size = bin_size
    nrows =  motifs_df.shape[0]

    bins = range(-upstream, (downstream), bin_size)
    bin_len = len(bins)
    motifs_concat_df = pd.concat([motifs_df]*bin_len, ignore_index="TRUE")
    motifs_sorted_df = motifs_concat_df.sort_values(["chrom","start","end"])

    ## Copy the bin list that is deep copy:
    bin_start_list = bins[:]
    bin_end_list = []
    for each in bin_start_list:
        bin_end_list.append(each+bin_size)

    bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
    bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
    bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

    ## Combine the motifs df and bin df by cbind or column wise:
    temp_motifs_df = pd.concat([motifs_sorted_df.reset_index(), bin_concat_df], axis = 1)
    temp_motifs_df["motifs_midpoint"] = (temp_motifs_df["start"] + temp_motifs_df["end"])/2
    temp_motifs_df["motifs_midpoint"] = temp_motifs_df["motifs_midpoint"].round().astype(int)
    final_motifs_df = temp_motifs_df.loc[:,["chrom", "motifs_midpoint", "bin_start", "bin_end", "start", "end", "motif_id", "tf_name"]]

    """ 
    chrom_start = tss_midpt + (bin_start); 
    chrom_end = tss_midpt + (bin_end) """
    final_motifs_df["chrom_start"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_start"]
    final_motifs_df["chrom_end"] = final_motifs_df["motifs_midpoint"] + final_motifs_df["bin_end"]

    select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end']
    final_motifs_df = final_motifs_df.loc[:,select_cols]

    ## Case handling with coordinates(mid point) smaller than the upstream range (-2000) indicated; else start coordinate < 0(-ve)
    final_motifs_df = final_motifs_df.loc[final_motifs_df["chrom_start"] > 0, :] 
    return(final_motifs_df)


def load_cpg_pybedtool_object(file_name_with_full_path):
    print "Processing Cpg bed file\n "
    cpg_df =  pd.read_csv(file_name_with_full_path, sep="\t", header=None)
    #cpg_df = cpg_df.iloc[:,0:6]
    Cpg_bed_file = pybedtools.BedTool.from_dataframe(cpg_df)
    return(Cpg_bed_file)


def generate_peaks_binned_perc_meth(peaks_final_bedfile, cpg_bedfile_list, **kwargs):
    # file_name =  kwargs["files_basename"]
    print "kwargs: ", kwargs
    peaks_bedfile = peaks_final_bedfile

    master_dict = {}
    for idx, each_file in enumerate(cpg_bedfile_list):
        cpg_bed_file = each_file
        cpg_bed_intersect = cpg_bed_file.intersect(peaks_bedfile, wa = True, wb = True)     
        cpg_bed_intersect = cpg_bedfile.intersect(peaks_bed_file, wa = True, wb = True)     
        # print(cpg_bed_intersect.head())  

        ## Working with the dataframes; reading the output file of pybedtool intersect:
        df_file = pd.read_csv(cpg_bed_intersect.fn, sep = "\t", header = None)
        df_ordered = df_file.iloc[:, 0:12]
        df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "wgbs_meth", "wgbs_unmeth", "chipbs_meth", "chipbs_unmeth", \
                            "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end"]
        df_ordered[["wgbs_meth", "wgbs_unmeth"]] = df_ordered[["wgbs_meth", "wgbs_unmeth"]].astype(int)
        df_ordered[["chipbs_meth", "chipbs_unmeth"]] = df_ordered[["chipbs_meth", "chipbs_unmeth"]].astype(int)
        df_ordered.to_csv(join(temp_output, sample + "_" + key + "_pyintersect.bed"), sep ="\t", header = True, index = False)
        df_ordered["count"] = 1
        df_grouped =  df_ordered.groupby(["bin_start", "bin_end"]).apply(lambda x : (x["chipbs_meth"].sum()/float(x["chipbs_meth"].sum() + x["chipbs_unmeth"].sum()),\
                                                                                                                 x["wgbs_meth"].sum()/float(x["wgbs_meth"].sum() + x["wgbs_unmeth"].sum()),
                                                                                                                 x["count"].sum()
                                                                                                                 ))   
        df_grouped_reset = df_grouped.reset_index(name="percent_info")
        df_grouped_reset[["chip_meth_perc","wgbs_meth_perc", "cpg_count_in_bin"]] = df_grouped_reset["percent_info"].apply(pd.Series)
        df_grouped_reset["cpg_count_in_bin"] = df_grouped_reset["cpg_count_in_bin"].astype(int)
        file_name = sample + "_" + key
        if file_name not in master_dict:
            master_dict[file_name] = df_grouped_reset
        print "Dimension of currently intersected peak file is", df_grouped_reset.shape

    ## Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
    combined_df = pd.concat(master_dict).reset_index().drop(["level_1"],1)
    combined_df.rename(columns={"level_0" : "annotation"}, inplace=True)
    combined_df["annotation"] = combined_df["annotation"].str.replace("_human|_FLAG", "")
    return(combined_df)


## For Single Cell-type : unbinned whole peak meth profile:
def main():
    ## Cell-type to analyse:
    #sample = "K562"
    #tf_list = ["ZBTB33", "YY1", "POLR2A", "TAF1", "SP1"]

    sample = "HepG2"
    tf_list = ["ZBTB33", "JUND", "GATA2", "POLR2A", "TAF1", "YY1"]

    ## Set output dir:
    output_dir = "/gpfs/gpfs1/home/schhetri/DNAme_paper/chipbs"

    ## motif peak file:
    fimo_file_path = "/gpfs/gpfs1/home/schhetri/DNAme_paper/unique_TF_files/" + \
                        sample + "/motif_calls/*/meme_fimo_dir/fimo.txt"

    ## HepG2 CpG file:
    cpg_file_path = "/gpfs/gpfs1/home/schhetri/DNAme_paper/chipbs/{}/merged_CpG_files/wgbs_chipbs_merged_CpG_files/" + \
                    "{}"  + "_" + "{}" + "_merged_CpG_context_deduplicated.bismark.cov.srt"


    ## Generate fimo motif file_dict:
    fimo_file_list = glob(join(fimo_file_path))
    fimo_filedict = {file.rsplit("/",3)[1].replace("_human", "") :file for file in fimo_file_list}
    chipbs_fimo_filedict = {tf:fimo_filedict.get(tf) for tf in tf_list}
    filt_fimo_filedict = {k:v for k,v in chipbs_fimo_filedict.iteritems() if v is not None}
    # chipbs_fimo_filedict.clear(); chipbs_fimo_filedict.update(filt_fimo_filedict) #retains ori dict

    ## Generate output dirs:
    dname_motif_output = join(output_dir, "motif_analysis", "motif_binned_methprofiles") 
    plot_output = join(dname_motif_output, sample + "_composite_plots")
    plot_data_output = join(dname_motif_output, sample + "_composite_plots_data")
    temp_output = join(dname_motif_output, "temp_cpg_motif")

    if not os.path.exists(dname_motif_output):
        os.makedirs(dname_motif_output)

    if not os.path.exists(plot_output):
        os.makedirs(plot_output) 

    if not os.path.exists(plot_data_output):
        os.makedirs(plot_data_output) 

    if not os.path.exists(temp_output):
        os.makedirs(temp_output) 

    ## Load CpG file only once to save computation time on each intersect
    concat_list = []
    for key, value in filt_fimo_filedict.iteritems():
        print("\nProcessing {} {}".format(sample, key))
        peaks_coord_df = final_fimo_motifs_model(filt_fimo_filedict.get(key), key, motif1=1, motif2=2)
        peaks_coord_df = peaks_coord_df.drop(["p-value"],1)
        peaks_binned_coord_df= generate_motifs_binned_coords(peaks_coord_df, 20, 20, 1)
        peaks_bed_file = pybedtools.BedTool.from_dataframe(peaks_binned_coord_df)
        cpg_bedfile = load_cpg_pybedtool_object(cpg_file_path.format(sample, sample, key))
        cpg_pybed_list = [cpg_bedfile]
        plot_data_set = generate_peaks_binned_perc_meth(peaks_bed_file, cpg_pybed_list)
        concat_list.append(plot_data_set)
        plot_data_set.to_csv(join(plot_data_output, key + "_binned_motif_methpercent_states.txt"), sep ="\t", header = True, index = False)

    final_meth_df = pd.concat(concat_list) #ignore_index=True
    final_meth_df.to_csv(join(dname_motif_output, sample + "_binned_motif_methpercent_states.txt"), sep ="\t", header = True, index = False)
    # final_meth_df_ = final_meth_df.pivot(index="sample", columns="bin_start", values="meth_percent")
    # final_meth_df_.to_csv(join(dname_motif_output, sample + "_motif_bin_methplot_heatmap_data.txt"), sep ="\t", header = True, index = True)

if __name__ == '__main__':
    main()

print("Time for analysis : {}".format(time.time()-start_time))
print("Task completed")


## Needed for multiple threads running with a race condition to create the dir:
# try:
#     os.makedirs(output_dir)
# except OSError as exc:
#     if exc.errno != errno.EEXIST:
#         raise
#     pass


