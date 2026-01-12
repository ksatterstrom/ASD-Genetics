# Do it with alternate rate files (these did not include OS variants; one was not depth filtered, the other was)

rates = read.table('ASD_bespoke_mutation_rates_draft_noOS_2024-07-30.txt',
                   header=T, sep='\t', stringsAsFactors = F) # 29,379
rates$chrom = gsub('chr', '', rates$chrom)

# Any duplicate genes?
nrow(rates[duplicated(rates$gene_id),]) # 0
nrow(rates[duplicated(rates$gene),])    # 20

# Define function

# Fixing genes where one gene ID maps to multiple ENSG IDs and vice versa
fix_mutation_rate_input <- function(rate_input) {
  rate_input$gene_id[rate_input$gene == 'HIST1H4F'] = 'ENSG00000274618' # (Note this is not a double here)
  rate_input$gene_id[rate_input$gene == 'ATXN7'] = 'ENSG00000163635'
  rate_input$gene_id[rate_input$gene == 'PRSS50'] = 'ENSG00000283706'
  rate_input$gene_id[rate_input$gene == 'CCDC39'] = 'ENSG00000284862'
  rate_input$gene_id[rate_input$gene == 'HSPA14'] = 'ENSG00000187522'
  rate_input$gene_id[rate_input$gene == 'IGF2'] = 'ENSG00000167244'
  rate_input$gene_id[rate_input$gene == 'PDE11A'] = 'ENSG00000128655'

  rate_input$gene[rate_input$gene_id == 'ENSG00000207258'] = 'RF00019_chr8'
  rate_input$gene[rate_input$gene_id == 'ENSG00000251792'] = 'RF00019_chr14'
  rate_input$gene[rate_input$gene_id == 'ENSG00000199668'] = 'RF00019_chr16'
  rate_input$gene[rate_input$gene_id == 'ENSG00000199444'] = 'RF00019_chr19'

  rate_input$gene[rate_input$gene == 'RF00019' & rate_input$chrom == '3'] = 'RF00019_chr3'
  rate_input$gene_id[rate_input$gene == 'RF00019_chr3'] = 'ENSG00000212392'

  rate_input$gene_id[rate_input$gene == 'ABCF2'] = 'ENSG00000285292'
  rate_input$gene_id[rate_input$gene == 'ELFN2'] = 'ENSG00000166897'
  rate_input$gene_id[rate_input$gene == 'GOLGA8M'] = 'ENSG00000188626'
  rate_input$gene_id[rate_input$gene == 'RGS5'] = 'ENSG00000143248'
  rate_input$gene_id[rate_input$gene == 'TMSB15B'] = 'ENSG00000158427'
  rate_input$gene_id[rate_input$gene == 'ZNF883'] = 'ENSG00000285447'

  rate_input$gene[rate_input$gene_id == 'ENSG00000232527'] = 'LSP1P5_chr1'

  rate_input$gene_id[rate_input$gene == 'HERC2P7'] = 'ENSG00000281909'
  rate_input$gene_id[rate_input$gene == 'LINC00484'] = 'ENSG00000235641'
  rate_input$gene_id[rate_input$gene == 'LINC01422'] = 'ENSG00000223704'
  rate_input$gene_id[rate_input$gene == 'LINC01505'] = 'ENSG00000234323'
  rate_input$gene_id[rate_input$gene == 'LINC02203'] = 'ENSG00000280709'
  rate_input$gene_id[rate_input$gene == 'LINC02528'] = 'ENSG00000229922'

  rate_input$gene[rate_input$gene_id == 'ENSG00000225255'] = 'LINC01297_chr22'
  rate_input$gene[rate_input$gene_id == 'ENSG00000274827'] = 'LINC01297_chr14'

  # Collect
  rate_input = group_by(rate_input, gene_id, gene, chrom)
  rate_input = summarise(rate_input, mu_snp_PTV = sum(mu_snp_PTV),
                           mu_snp_Mis2 = sum(mu_snp_Mis2),
                           mu_snp_Mis1 = sum(mu_snp_Mis1),
                           mu_snp_Mis0 = sum(mu_snp_Mis0),
                           mu_snp_Syn = sum(mu_snp_Syn),
                           mu_snp_MisB = sum(mu_snp_MisB),
                           mu_snp_MisA = sum(mu_snp_MisA),
                           mu_snp_MisOther = sum(mu_snp_MisOther) )

  return(data.frame(rate_input))
}

# Apply function
rates = fix_mutation_rate_input(rates)

# Which genes are relevant to this analysis?

# Read all genes involved here:
gene_table = read.table('ASD_VEP_gene_list_2024-07-12.txt',
                        header=T, sep='\t', stringsAsFactors = F) # 21,777

nrow(subset(rates, gene %in% gene_table$Gene))       # 21,138
nrow(subset(rates, gene_id %in% gene_table$Gene_ID)) # 21,138

nrow(subset(gene_table, Gene %in% rates$gene))       # 21,138
nrow(subset(gene_table, Gene_ID %in% rates$gene_id)) # 21,138

gene_table = merge(gene_table, rates, by.x='Gene_ID', by.y='gene_id', all.x=T, all.y=F) # 21,777

gene_table$mu_snp_MisAll = gene_table$mu_snp_MisB + gene_table$mu_snp_MisA + gene_table$mu_snp_MisOther

# This comes from autosomal work (Fu et al., 2022)
gene_table$mu_total_PTV = 1.6 * gene_table$mu_snp_PTV

gene_table = gene_table[c('Gene', 'Gene_ID', 'Chrom', 'LOEUF', 'LOEUF_bin', 'mu_total_PTV', 
                            'mu_snp_MisB', 'mu_snp_MisA', 'mu_snp_MisAll',  'mu_snp_Syn')]

colnames(gene_table) = c('Gene', 'Gene_ID', 'Chrom', 'LOEUF', 'LOEUF_bin', 'mu_PTV',
                          'mu_MisB', 'mu_MisA', 'mu_MisAll', 'mu_Syn')

gene_table = gene_table[order(gene_table$Gene),]

write.table(gene_table, 'ASD_gene_table_w_bespoke_mutation_rates_noOS_yes_covg_filter_2024-07-30.txt',
            sep='\t', row.names=F, col.names=T, quote=F)