#! /bin/bash
# Using pdfjam to combine figures into single doc
# working command: pdfjam s-alpinus-boxplots.pdf salmo-salar-boxplots.pdf --nup 1x2 --outfile combo.pdf
pdfjam ./outputs/108/t-thymallus-boxplots.pdf ./outputs/108/salmo-salar-boxplots.pdf \
./outputs/108/s-alpinus-boxplots.pdf ./outputs/108/o-mykiss-boxplots.pdf \
./outputs/108/o-tshaw-boxplots.pdf ./outputs/108/o-kisutch-boxplots.pdf \
--nup 1x6 --outfile ./outputs/109/Fig4.pdf