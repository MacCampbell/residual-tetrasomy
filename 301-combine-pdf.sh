#! /bin/bash
# Using pdfjam to combine figures into single doc
# working command: pdfjam s-alpinus-boxplots.pdf salmo-salar-boxplots.pdf --nup 1x2 --outfile combo.pdf
pdfjam ./outputs/201/t-thymallus-boxplots.pdf ./outputs/201/salmo-salar-boxplots.pdf \
./outputs/201/s-alpinus-boxplots.pdf ./outputs/201/o-mykiss-boxplots.pdf \
./outputs/201/o-tshaw-boxplots.pdf  \
--nup 1x5 --outfile ./outputs/301/Fig4.pdf

pdfjam ./outputs/201/t-thymallus-predictions.pdf ./outputs/201/salmo-salar-predictions.pdf \
./outputs/201/s-alpinus-predictions.pdf ./outputs/201/o-mykiss-predictions.pdf \
./outputs/201/o-tshaw-predictions.pdf \
--nup 1x5 --outfile ./outputs/301/S8_support_for_assignment.pdf