#!/bin/bash

scp -pr annotations/ salias@etna:/home/salias/TFM/
scp -pr coex-analysis/ salias@etna:/home/salias/TFM/
scp COTAN/results/COTAN-cluster-script.R salias@etna:/home/salias/TFM/COTAN/results
scp COTAN/results/coex-to-tsv.R salias@etna:/home/salias/TFM/COTAN/results
