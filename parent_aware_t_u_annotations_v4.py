# Parent-aware TDT annotations
#
# This draws heavily from Jack Kosmicki's function, 
# https://hail.is/docs/0.2/methods/genetics.html#hail.methods.transmission_disequilibrium_test,
# and, when summed, should give the same output values
#
# Requires trio dataset and, like Hail's TDT function, does not cover the Y
#
# Also note:
# 1) Uses sexes from ped and assumes fathers are male and mothers are female
# 2) Does not allow fathers to be het when they should be hemizygous
# 3) To match Jack's function, requires fathers to have a genotype even when considering regions where 
#    proband is hemizygous
def parent_aware_t_u_annotations_v4(td):

    # First decide copy state
    td = td.annotate_entries(autosomal_copy_state = ( hl.case()
            .when(td.locus.in_autosome() | td.locus.in_x_par() | (td.is_female == True), True)
            .when(td.locus.in_x_nonpar() & (td.is_female == False), False)
            .default(hl.missing('bool')) ) )
    # Note: the above uses the "is_female" from the ped and not from the dataset itself
    
    # Now annotate t & u values
    td = td.annotate_entries(
        
        t_from_dad = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
            ( (td.proband_entry.GT.is_het() & td.mother_entry.GT.is_hom_ref()) | 
              (td.proband_entry.GT.is_hom_var() & 
               (td.mother_entry.GT.is_het() | td.mother_entry.GT.is_hom_var()) )) ), 1, 0),
        
        t_from_mom = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.mother_entry.GT.is_het() &
            ( (td.proband_entry.GT.is_het() & td.father_entry.GT.is_hom_ref()) | 
              (td.proband_entry.GT.is_hom_var() & 
               ((td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False)) | td.father_entry.GT.is_hom_var()) )) ) |
           ( (td.autosomal_copy_state == False) & td.mother_entry.GT.is_het() &
               td.proband_entry.GT.is_hom_var() & 
               (td.father_entry.GT.is_hom_ref() | td.father_entry.GT.is_hom_var()) ), 1, 0),
        # I could consider removing any reference at all to father's genotype in this last line    
        
        u_from_dad = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
            ( (td.proband_entry.GT.is_het() & td.mother_entry.GT.is_hom_var()) | 
              (td.proband_entry.GT.is_hom_ref() & 
               (td.mother_entry.GT.is_het() | td.mother_entry.GT.is_hom_ref()) )) ), 1, 0),
        
        u_from_mom = hl.if_else(
           ( (td.autosomal_copy_state == True) & td.mother_entry.GT.is_het() &
            ( (td.proband_entry.GT.is_het() & td.father_entry.GT.is_hom_var()) | 
              (td.proband_entry.GT.is_hom_ref() & 
               ((td.father_entry.GT.is_het() & (td.locus.in_x_nonpar() == False)) | td.father_entry.GT.is_hom_ref()) )) ) |
           ( (td.autosomal_copy_state == False) & td.mother_entry.GT.is_het() &
               td.proband_entry.GT.is_hom_ref() & 
               (td.father_entry.GT.is_hom_ref() | td.father_entry.GT.is_hom_var()) ), 1, 0),   
        # Again, could consider removing any reference at all to father's genotype in this last line
        
        t_indeterminate = hl.if_else( 
            (td.autosomal_copy_state == True) & td.proband_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
             td.father_entry.GT.is_het() & td.mother_entry.GT.is_het(), 1, 0),

        u_indeterminate = hl.if_else( 
            (td.autosomal_copy_state == True) & td.proband_entry.GT.is_het() & (td.locus.in_x_nonpar() == False) & 
             td.father_entry.GT.is_het() & td.mother_entry.GT.is_het(), 1, 0)        
    )
        
    return (td)
