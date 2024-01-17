#!/bin/bash

echo $1

rm -rf $1/fpr \
       $1/filters \
       $1/howde_index \
       $1/merge_infos \
       $1/counts \
       $1/histograms \
       $1/partition_infos \
       $1/superkmers \
       $1/minimizers
