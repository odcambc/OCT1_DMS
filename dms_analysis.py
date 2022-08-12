import process_oligo_list
import process_variants
import csv
import pandas as pd
from Bio import SeqIO

import yaml


#TODO: check for 4x deletions - analyze those as well!

order = ["A", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
         "D_1", "D_2", "D_3",
         "I_1", "I_2", "I_3", "X"]

# Oligo processing files

designed_variant_header = ["count", "pos", "chunk_pos", "chunk", "mutation_type", "name", "codon", "mutation", "length",
                           "hgvs"]

# The reference fasta files include pre-ORF nucleotides. This offset value accounts for that.

def process_experiment(experiment_yaml, regenerate_variants = False):

    with open(experiment_yaml, 'r') as file:
        experiment_definition = yaml.safe_load(file)

        ref_file = experiment_definition['ref_file']
        ref_name = experiment_definition['ref_name']
        oligo_file = experiment_definition['oligo_file']
        offset = experiment_definition['offset']
        variants_file = experiment_definition['variants_file']
        chunk_starts = experiment_definition['chunk_starts']
        chunk_ends = experiment_definition['chunk_ends']
        gatk_output_dir = experiment_definition['gatk_output_dir']
        files = experiment_definition['files']
        experiments = experiment_definition['experiments']
        experiment_name = experiment_definition['experiment_name']

        chunks = list(zip(chunk_starts, chunk_ends))

    with open(ref_file) as f:
        ref_dict = SeqIO.to_dict(SeqIO.parse(ref_file, "fasta"))

    if regenerate_variants:
        variant_list = process_oligo_list.designed_variants(oligo_file, str(ref_dict[ref_name].seq),
                                                                offset,
                                                                chunks)
        process_oligo_list.write_designed_csv(variants_file, designed_variant_header, variant_list)

    with open(variants_file, 'r') as f:
        variants_reader = csv.reader(f, delimiter=',')
        designed_variants_list = list(variants_reader)

    designed_df = pd.DataFrame.from_records(designed_variants_list[1:],
                                                  columns=designed_variants_list[0]).convert_dtypes()
    designed_df["pos"] = pd.to_numeric(designed_df["pos"])
    designed_df["count"] = pd.to_numeric(designed_df["count"])
    designed_df["chunk_pos"] = pd.to_numeric(designed_df["chunk_pos"])
    designed_df["chunk"] = pd.to_numeric(designed_df["chunk"])

    for file, experiment in zip(files, experiments):
        csv_file = process_variants.read_gatk_csv(gatk_output_dir + experiment_name + "/" + file)
        df, other = process_variants.process_variants_file(csv_file, designed_df)

        # Write an Enrich2-readable output
        process_variants.write_enrich_df("data/" + experiment_name + "/" + experiment + ".tsv", df)

        # Write the processed file to tsv
        df.to_csv("data/counts/" + experiment_name + "/" + experiment + ".csv")

        # Write the rejected variants as well
        with open("data/" + experiment_name + "/rejected/rejected_" + experiment + ".csv", 'w') as f:
            csvwriter = csv.writer(f)
            csvwriter.writerows(other)

    # Remove any variants with no observations before processing with Enrich2.
    if remove_zeros:
        enrich_file_list = ["data/" + experiment_name + "/" + s + ".tsv" for s in experiments]
        zeros = process_variants.remove_zeros_enrich(enrich_file_list)

        with open("data/" + experiment_name + "/rejected/" + experiment_name + "_unobserved_variants.csv", 'w+') as f:
            for variant in zeros:
                f.write("%s\n" % variant)

if __name__ == '__main__':

    TRPV1 = False
    Kir = True
    VatD = False
    OCT1 = False
    test = False
    pacbio = False
    remove_zeros = True
    kir21_AEZ = False
    baseline = False
    OCT1_full = False
    MOR = False
    TRPV1_MiSeq = False

    #process_experiment('config_files/kir21_full.yaml')
    #process_experiment('config_files/mor_library.yaml', regenerate_variants = False)
    #process_experiment('config_files/trpv1_novaseq_qc.yaml', regenerate_variants = False)
    #process_experiment('config_files/vatd.yaml', regenerate_variants = False)
    process_experiment('config_files/met.yaml', regenerate_variants = False)