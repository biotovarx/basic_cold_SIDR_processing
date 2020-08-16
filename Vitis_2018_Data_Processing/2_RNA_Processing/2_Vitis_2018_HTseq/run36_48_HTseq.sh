#!/bin/bash



# Run 36 hour bash

bash ~/data_processing/htseq_counting/Vit_36hr/runHTseq_BatchInput_A.sh ~/data_processing/star_alignment/Vitis_2018_STAR/Vitis_TR_QT_QF_QF_2018_STAR_Genome_2018-08-03_090538/Vitis_36hpc_*

# Run 48 hour bash

bash ~/data_processing/htseq_counting/Vit_48hr/runHTseq_BatchInput_A.sh ~/data_processing/star_alignment/Vitis_2018_STAR/Vitis_TR_QT_QF_QF_2018_STAR_Genome_2018-08-03_090538/Vitis_48hpc_*
