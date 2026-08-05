#!/usr/bin/env bash
set -euo pipefail

Rscript -e "rmarkdown::render('index.Rmd')"
