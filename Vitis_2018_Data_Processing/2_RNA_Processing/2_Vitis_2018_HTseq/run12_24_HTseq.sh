#!/bin/bash



# Run 12 hour bash

bash ~/data_processing/htseq_counting/Vit_12hr/runHTseq_BatchInput_A.sh ~/data_processing/star_alignment/Vitis_2018_STAR/Vitis_TR_QT_QF_QF_2018_STAR_Genome_2018-08-03_090538/Vitis_12hpc_*

# Run 24 hour bash

bash ~/data_processing/htseq_counting/Vit_24hr/runHTseq_BatchInput_A.sh ~/data_processing/star_alignment/Vitis_2018_STAR/Vitis_TR_QT_QF_QF_2018_STAR_Genome_2018-08-03_090538/Vitis_24hpc_*
