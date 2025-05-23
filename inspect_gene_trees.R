library(ape)
source("funcs.R")
args = commandArgs(trailingOnly=TRUE)
if(length(args)<5){
  msg = c("Usage: Rscript gene_trees.R <BRLEN_THRESH> <TAXA_THRESH> <TREEFILE> <GENE_NAMES_FILE> <OUTPUT_DIR_NAME>\n",
          "==============================================================================================================\n",
          "\tBRLEN_THRESH: maximum relative branch length of the longest branch (0-1)\n",
          "\tTAXA_THRESH: minimum proportion of taxa needs to be retained by the 'saved' gene compared to the original gene tree (0-1)\n",
          "\tTREEFILE: name to file that holds all the gene trees (one tree per line)\n",
          "\tGENE_NAMES_FILE: file with one gene name per line\n",
          "\tOUTPUT_DIR_NAME: name of output dir to save results from a single run")
  msg = paste0(msg, collapse="")
  stop(msg)
}
print(args)
threshold = args[1]
tax_threshold = args[2]
trees_file = args[3]
gene_names_file = args[4]
output_dir=args[5]

if(!(dir.exists(output_dir))){
  dir.create(output_dir)
}
### Reading in the tree and the data -----
potential_outgroups = c(
  "GCA_958295455",
  "GCA_027580185",
  "GCA_022674325",
  "GCA_033807575"
)
gene_trees = read.tree(trees_file)
gene_trees = lapply(gene_trees, function(tree){
  return(tryCatch({ey.reroot_phylo(tree, potential_outgroups)},error = function(e){ 
      print(e)
      return(tree)
    }
  ) 
  )
})
gene_names = scan(gene_names_file, what=character())
gene_names = sapply(strsplit(gene_names, "/"), function(x){return(x[2])})
#setwd(output_dir)

max_n_tax = max(sapply(gene_trees, function(x){length(x$tip.label)}))

### Calculating the relative length of the longest branch ----
brlens = lapply(gene_trees,function(X){X$edge.len})
prop_max = sapply(brlens, function(X){max(X)/sum(X)})
names(prop_max) = gene_names
names(gene_trees) = gene_names

prop_tax = sapply(gene_trees, function(x){length(x$tip.label)/max_n_tax})
names(prop_tax) = gene_names
### Filtering and processing trees according to defined thresholds ----

long_branch_genes = prop_max[prop_max > threshold]
long_branch_enough_taxa_genes = names(long_branch_genes)

#### Running tree splitting function on suitable genes -----

saved_trees = sapply(long_branch_enough_taxa_genes, function(gene){
    this_tree = gene_trees[[gene]]
    this_tree_saved = split_trees(this_tree, retain_threshold = tax_threshold, prop_max_threshold = threshold)
    if(!is.null(this_tree_saved)){
        if(Ntip(this_tree_saved)/Ntip(this_tree) < tax_threshold){
            this_tree_saved = NULL
        }
    }
    return(this_tree_saved)
})
pdf(paste(output_dir,"saved_genes.pdf", sep="/"), 20,10)
if(length(saved_trees) !=0) {
    sum(sapply(saved_trees, is.null))
    par(mfrow=c(1,2))
    sapply(names(saved_trees), function(gene){
        orig_tree = gene_trees[[gene]]
        saved_tree = saved_trees[[gene]]
        plot(orig_tree, cex=0.5, main=gene)
        if(is.null(saved_tree)){
            plot.new()
        }else{
            plot(saved_tree, cex=0.5, main=paste(gene,"_saved", sep=""))
        }
        return(NULL)
    })
}
dev.off()


if(length(saved_trees) !=0) {
    saved_genes = names(saved_trees)[sapply(saved_trees, function(x){!is.null(x)})]
}

outlier_genes = names(prop_max[prop_max > threshold | prop_tax < tax_threshold])

if(length(saved_trees) !=0) {
    outlier_genes = setdiff(outlier_genes, saved_genes)
}

pdf(paste(output_dir,"outlier_genes.pdf", sep="/"),30,30)
par(mfrow=c(5,5))
sapply(outlier_genes, function(X){
    plot(gene_trees[[X]], main=X)
})
dev.off()
write(outlier_genes, file = paste(output_dir,"outlier_genes.txt", sep="/"),sep='\n')

#### Write out kept taxa for saved genes ------
    out.dir = paste(output_dir,"saved_genes_kept_taxa", sep="/")
    if (!file.exists(out.dir)){
        dir.create(out.dir)
    }

if(length(saved_trees) !=0) {
    sapply(saved_genes, function(gene){
        nm = paste(gene, "_kept_taxa.txt", sep="")
        nm = file.path(out.dir, nm)
        taxa = saved_trees[[gene]]$tip.label
        write(taxa, file=nm, sep='\n')
    })
}


### Write out non-outlier genes
outlier_or_saved = c(saved_genes, outlier_genes)
OK_genes = setdiff(gene_names, outlier_or_saved)

pdf(paste(output_dir,"OK_genes.pdf", sep="/"),30,30)
par(mfrow=c(5,5))
sapply(OK_genes, function(X){
    plot(gene_trees[[X]], main=X)
})
dev.off()
