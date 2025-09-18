# These are the functions I used to do the c-alpha tests comparing the distribution of variants across genes in one sample group to that of another sample group
# 
# These functions are built using the c-alpha functions in the R "Assotester" package as a starting point
# Those functions condition p0 on sample number; here I condition on variant number
#
# When I call this, it looks like:
# variant_list = subset(ptv_mis2_variants, (Role == 'Proband'))
# variant_list$Pheno = variant_list$DDID
#
# set.seed(42)
# calpha_PTV = calpha_function_newer(variant_list, perms, singletons = 'keep')

calpha_method_newer <- function(variant_list_with_phenotypes, p0, singletons) {
  # ASSUMES ALL VARIANTS ARE HETS AND THERE IS A BINARY "PHENO" COLUMN
  
  # Variant counts by gene and phenotype
  gene_table = group_by(variant_list_with_phenotypes, Gene)
  gene_table = summarise(gene_table, case_alleles = sum(Pheno), total_alleles = sum(nNonRef))
  
  if(singletons == 'group') {
    singles = subset(gene_table, total_alleles == 1)
    gene_table = subset(gene_table, total_alleles > 1)
    
    gene_table = rbind(gene_table, data.frame(Gene = 'SINGLETONS', 
                                              case_alleles = sum(singles$case_alleles), 
                                              total_alleles = sum(singles$total_alleles)))
  }
  
  # Test statistic
  Talpha = sum((gene_table$case_alleles - p0*gene_table$total_alleles)^2 - (gene_table$total_alleles * p0 * (1-p0)))
  
  return(Talpha)
}

calpha_function_newer <- function(variant_list_with_phenotypes, perm, singletons) {
  # ASSUMES ALL VARIANTS ARE HETS AND THERE IS A BINARY "PHENO" COLUMN
  
  start.time = Sys.time()
  
  if (!singletons %in% c('keep', 'group', 'drop')) { print ('No special handling for singletons by default') }
  
  # Variant counts by gene
  gene_table = table(variant_list_with_phenotypes$Gene)
  
  # Group singletons?
  if(singletons == 'group') {
    singles = gene_table[gene_table == 1]
    gene_table = gene_table[gene_table > 1]
    
    gene_table = c(gene_table, SINGLETONS = sum(singles))
  }
  
  # Drop singletons?
  if(singletons == 'drop') {
    gene_table = gene_table[gene_table > 1]
    
    # Adjust p0
    variant_list_with_phenotypes = subset(variant_list_with_phenotypes, Gene %in% rownames(gene_table))
  }
  
  # Get number of genes
  num_genes = length(gene_table)
  
  # Calculate p0...
  # How many alleles from cases
  nAff = sum(variant_list_with_phenotypes$Pheno)
  # How many total alleles
  nTot = nrow(variant_list_with_phenotypes)
  # Proportion of alleles from cases
  p0 = nAff / nTot  
  
  # Test statistic 
  calpha.stat = calpha_method_newer(variant_list_with_phenotypes, p0, singletons)
  
  # Variance of Talpha
  Valpha = 0
  for (i in 1:num_genes) {
    for (u in 0:gene_table[i]) {
      Valpha = Valpha + (((u - gene_table[i]*p0)^2 - gene_table[i]*p0*(1-p0))^2)*dbinom(u, gene_table[i], p0)
    }
  }
  names(Valpha) = NULL
  
  # Asymptotic p-value
  if (Valpha==0) { 
    asym.pval = 1 
  } else {
    asym.pval = 1 - pchisq(calpha.stat^2 / Valpha, df=1)
  }
  
  # Permutations
  perm.pval = NA
  
  x.perm = rep(0, perm)
  for (i in 1:perm) {
    permuted_variant_list = variant_list_with_phenotypes
    permuted_variant_list$Pheno = sample(permuted_variant_list$Pheno)
    x.perm[i] = calpha_method_newer(permuted_variant_list, p0, singletons) 
  }
  # p-value 
  perm.pval = sum(x.perm^2 > calpha.stat^2) / perm
  
  ## results
  name = "CALPHA: c-alpha Test"
  arg.spec = c(nAff, nTot-nAff, num_genes, perm)
  names(arg.spec) = c("case.alleles", "control.alleles", "genes", "n.perms")
  res = list(calpha.stat = calpha.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name)
  class(res) = "assoctest"	
  
  end.time = Sys.time()
  print(end.time - start.time)
  
  print(res)
  
  return(res)
}
