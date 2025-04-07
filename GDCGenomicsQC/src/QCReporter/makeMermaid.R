# formatted strings will be easier to update the numbers with
# Wherever you need to put a number put %s
# Then the following arguments replace the %s in order
# First the header

df <-
    data.frame(
        inSubjects = 1,
        inSNPs = 2,
        filteredSubjects = 3,
        filteredSNPs = 4,
        joinedSubjects = 5,
        joinedSNPs = 6, mafSNPs = 7,
        snpsMissing = 8,
        subjMissing = 9
    )


sprintf(
    "flowchart TB
    input(%s Subjects %s SNPs)
    filt(\"MAF > 0.05, (%s)  SNP Missingness < 0.05 (%s)  Subj Missing >0.05, (%s)\")
    joinedData(%s Subjects  %s SNPs)
    input --> filt
    filt --> joinedData",
    df$inSubjects,
    df$inSNPs,
    df$mafSNPs,
    df$snpsMissing,
    df$subjMissing,
    df$joinedSubjects,
    df$joinedSNPs
    
) |>
    writeLines("filtSteps.mmd")
