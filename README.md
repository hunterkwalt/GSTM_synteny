**get_gstm_genes.sh:** Bash script to parse GFFs and return flanking genes. This should be run in a directory containing folders named by your taxa of interest. Inside each folder should be the organisms gff annotation. The first 35 lines of the script will work for any RefSeq GFF. The remaining lines do the same thing, but I had to do a little customization for the Neotoma and Phyllotis annotation files.

**GSTm_synteny.R** R script to visualize results from get_gstm_genes.sh. This can be run with **gene_flank_locs_final.tsv** as input
