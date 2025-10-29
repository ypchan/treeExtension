# treeExtension

Tiny R utilities for working with phylogenetic trees and IQ-TREE outputs.

## Features

- Helpers to convert and annotate IQ-TREE files
- Utilities to map BLASTN results onto IQ-TREE trees

## Installation

```r
# if needed:
install.packages("remotes")

# install the package from GitHub
remotes::install_github("ypchan/treeExtension")
```

## Quick start

```r
library(treeExtension)
# see package help
help(package = "treeExtension")
```

## Development

Run these steps from the package root (recommended: open project in RStudio):

```r
# set working directory to the package root when using this.path
library(this.path)
setwd(this.dir())

# develop
library(devtools)
devtools::document()    # regenerate documentation
devtools::load_all()    # load package without reinstalling

# view function help
?make_binary_iqtree
?map_blastn_2_iqtree
```

## Key functions

- make_binary_iqtree: create binary annotations for IQ-TREE tree files.
- map_blastn_2_iqtree: map BLASTN results onto IQ-TREE trees for visualization/annotation.

Contributions and issues: open a GitHub issue or pull request in the repository.
