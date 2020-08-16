#!/bin/bash



# Run 12 hour bash

bash ~/data_processing/htseq_counting/Athaliana_2018_HTseq/AT_12hr/runHTseq_BatchInput_A.sh ~/data_processing/star_alignment/Arabidopsis_2018_STAR/Arabidopsis_AF_TR_QT_QF_QF_2018_STAR_Genome_2018-07-21_145716/At_*_12hpc*

# Run 24 hour bash

bash ~/data_processing/htseq_counting/Athaliana_2018_HTseq/AT_24hr/runHTseq_BatchInput_A.sh ~/data_processing/star_alignment/Arabidopsis_2018_STAR/Arabidopsis_AF_TR_QT_QF_QF_2018_STAR_Genome_2018-07-21_145716/At_*_24hpc*
