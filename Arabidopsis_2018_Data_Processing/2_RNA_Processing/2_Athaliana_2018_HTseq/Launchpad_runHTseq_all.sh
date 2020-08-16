#!/bin/bash


source activate ~/conda_temp_env/


nohup bash run12_24_HTseq.sh &


nohup bash run36_48_HTseq.sh &


nohup bash runUTC_HTseq.sh &


source deactivate
