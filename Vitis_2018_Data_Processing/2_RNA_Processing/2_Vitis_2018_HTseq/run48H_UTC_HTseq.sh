#!/bin/bash



# Run 48 hold bash

bash ~/data_processing/htseq_counting/Vit_48hold/runHTseq_BatchInput_A.sh ~/data_processing/star_alignment/Vitis_2018_STAR/Vitis_TR_QT_QF_QF_2018_STAR_Genome_2018-08-03_090538/Vitis_48hold_*

# Run UTC bash

bash ~/data_processing/htseq_counting/Vit_UTC/runHTseq_BatchInput_A.sh ~/data_processing/star_alignment/Vitis_2018_STAR/Vitis_TR_QT_QF_QF_2018_STAR_Genome_2018-08-03_090538/Vitis_UTC_*
