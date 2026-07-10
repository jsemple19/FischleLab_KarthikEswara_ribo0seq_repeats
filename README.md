# FischleLab_KarthikEswara_ribo0seq_repeats

## Preparing gtf

Nonredunant genome annotations for repeats were obtained fom [Dfam database version 3.5](https://www.dfam.org/releases/Dfam_3.5/annotations/ce10/) and lifted over to the ce11 version of the genome. Repeat loci were given unique names by appending "_rpt#" to the repeat family name, and a gtf was constructed with the same annotation given gene, transcript and exon types. This was then concatenated with canonical gene annotation gtf downloaded from [Wormbase (WS298 version)](https://downloads.wormbase.org/releases/WS298/species/c_elegans/PRJNA13758/). This was performed with the `getRepeatData.R` script.

Repeat classification was downloaded with [FamDB tools](https://github.com/Dfam-consortium/FamDB) for Dfam version 4.0, and then merged with Dfam 3.5 annotations using the `getRepeatAnnotation.R`script.



