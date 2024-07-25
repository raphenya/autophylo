import argparse, os, sys
from autophylo import tree_builder

def main():
    cpu_count=os.cpu_count()
    parser = argparse.ArgumentParser(prog="autophylo",description="Program to automatically build a phylogentic tree given a protein sequence(s)")
    parser.add_argument('-i','--input_sequence', dest="input_sequence", required=True, help='input file must be in either FASTA (protein)')
    parser.add_argument('-m','--multi_sequence', dest="multi_sequence", required=False, help='input file with multiple FASTA (protein)')
    parser.add_argument('-o','--output_file', dest="output_file", required=True, help="output folder and base filename")
    parser.add_argument('-c','--config', dest="config_file", required=True, help="path to config file")
    parser.add_argument('-s','--sample_size', dest="sample_size", required=False, default=250, help="sample size (default=250)")
    parser.add_argument('-n','--threads', dest="threads", type=int,default=cpu_count, help="number of threads (CPUs) to use (default={})".format(cpu_count))
    # clustering
    parser.add_argument('--skip_usearch', dest="skip_usearch", action="store_true", help="skip clustering usearch (default=False)")
    parser.add_argument('--usearch_cluster_size', dest="usearch_cluster_size", required=False, default=4, help="usearch cluster size (default=4)")

    # use cached files
    parser.add_argument('--skip_taxids', dest="skip_taxids", action="store_true", help="skip taxids download (default=False)")
    parser.add_argument('--skip_blast', dest="skip_blast", action="store_true", help="skip blast (default=False)")
    # TODO: to be removed
    parser.add_argument('--skip_sampling', dest="skip_sampling", action="store_true", help="skip sampling (default=False)")

    # filters
    parser.add_argument('--specific_taxa', dest="specific_taxa",type=str.lower, default="",
                        help = "specify taxonomy (default = 'actinobacteriota','bacteroidota','desulfobacterota','firmicutes','proteobacteria','synergistota','verrucomicrobiota','fusobacteria')")
    parser.add_argument('--minimum_length', dest="minimum_length", action="store_true", help="minimum sampled sequence length (default=False)")
    parser.add_argument('--maximum_length', dest="maximum_length", action="store_true", help="maximum sampled sequence length (default=False)")
    parser.add_argument('--percent_positive_scoring', dest="percent_positive_scoring", required=False, default=50.0, help="percent positive scoring (default=50.0)")
    parser.add_argument('--percent_identity', dest="percent_identity", required=False, default=50.0, help="percent positive scoring (default=50.0)")
    parser.add_argument('--usearch_percent_identity', dest="usearch_percent_identity", required=False, default=0.75, help="usearch_percent_identity (default=0.75)")
    parser.add_argument('--remove_similar_sequences', dest="remove_similar_sequences", action="store_true", help="remove similar sequences with 2 or 3 amino acid difference (default=False)")

    # ML options
    parser.add_argument('--skip_maximum_likelihood_tree', dest="skip_maximum_likelihood_tree", action="store_true", help="skip maximum likelihood tree generation (default=False)")
    parser.add_argument('--RAxML', dest="RAxML",type=str.upper,choices = ['PTHREADS', 'HYBRID'],default="PTHREADS",help = "specify flavour of RAxML to run (default = PTHREADS)")
    # specify model to use for ML trees
    parser.add_argument('--model', dest="model",default="PROTCATIJTTF",help = "specify Maximum Likelihood model to use (default = 'PROTCATIJTTF')")

    parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
    parser.add_argument('--clean', dest="clean", action="store_true", help="removes temporary files")

    if len(sys.argv) == 1:
        sys.stderr.write("No arguments provided, printing help menu ...\n")
        parser.print_help()
        sys.exit(1)
        
    args=parser.parse_args()
    tree_builder.main(args)

if __name__=="__main__":
    main()