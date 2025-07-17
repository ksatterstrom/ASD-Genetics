# Define sex-aware variant call rate calculation with Hardy-Weinberg p value
def sex_aware_variant_annotations_with_pHWE(mt):
    num_males = mt.aggregate_cols(hl.agg.count_where(mt.is_female == False))
    num_females = mt.aggregate_cols(hl.agg.count_where(mt.is_female == True))
    
    mt = mt.annotate_rows(
        male_hets = hl.agg.count_where(mt.GT.is_het() & (mt.is_female == False)),
        male_homvars = hl.agg.count_where(mt.GT.is_hom_var() & (mt.is_female == False)),
        male_calls = hl.agg.count_where(hl.is_defined(mt.GT) & (mt.is_female == False)),
        female_hets = hl.agg.count_where(mt.GT.is_het() & (mt.is_female == True)),
        female_homvars = hl.agg.count_where(mt.GT.is_hom_var() & (mt.is_female == True)),
        female_calls = hl.agg.count_where(hl.is_defined(mt.GT) & (mt.is_female == True))
    )
    
    mt = mt.annotate_rows(
        call_rate = ( hl.case()
            .when(mt.locus.in_y_nonpar(), (mt.male_calls / num_males))
            .when(mt.locus.in_x_nonpar(), (mt.male_calls + 2*mt.female_calls) / (num_males + 2*num_females))
            .default((mt.male_calls + mt.female_calls) / (num_males + num_females)) ),
        AC = ( hl.case()
            .when(mt.locus.in_y_nonpar(), mt.male_homvars)
            .when(mt.locus.in_x_nonpar(), mt.male_homvars + mt.female_hets + 2*mt.female_homvars)
            .default(mt.male_hets + 2*mt.male_homvars + mt.female_hets + 2*mt.female_homvars) ),
        AN = ( hl.case()
            .when(mt.locus.in_y_nonpar(), mt.male_calls)
            .when(mt.locus.in_x_nonpar(), mt.male_calls + 2*mt.female_calls)
            .default(2*mt.male_calls + 2*mt.female_calls) ),
        pHWE = ( hl.case() 
            .when(mt.locus.in_y_nonpar() | mt.locus.in_mito(), 1.0)
            .when(mt.locus.in_x_nonpar(), hl.hardy_weinberg_test(
                hl.int32(mt.female_calls - mt.female_hets - mt.female_homvars), 
                hl.int32(mt.female_hets), 
                hl.int32(mt.female_homvars)).p_value)
            .default(hl.hardy_weinberg_test(
                hl.int32(mt.male_calls+mt.female_calls - mt.male_hets-mt.female_hets - mt.male_homvars-mt.female_homvars), 
                hl.int32(mt.male_hets + mt.female_hets), 
                hl.int32(mt.male_homvars + mt.female_homvars)).p_value) )
    )
    
    return (mt)
