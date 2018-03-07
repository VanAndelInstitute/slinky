# Abstract

The L1000 data set from the NIH LINCS program holds the promise to deconvolute a wide range of biological questions in transcriptional 
space.  However, using this large and decentralized data set presents its own challenge.  The `slinky` package was created to to simplify filtering and accessing this data in an efficient analysis pipeline.

# Installation

We are in the process of preparing to submit this package to the BioConductor project.  In the mean time, you can install it and its dependencies as follows:

```
devtools::install_git("https://github.com/vanandelinstitute/slinky")
```

# Prerequisites

## Memory

Slinky uses the rhdf5 package to access slices of the L1000 data from disk.  As a result, you only need enough memory to hold the data slices you are working with, and possibly the metadata as well.  (The entire metadata data frame for all 1.3 million instances requires a modest 220MB of RAM). The development machine on which the slinky package was primarily developed has 16GB of RAM and can load and manipulate tens of thousands of instances from the L1000 dataset without difficulty.  Of course, YMMV depending on the specifications of your computer.

## Data

The examples in the package vignette can be completed with the demonstration `gctx` 
and `info` files that are installed along with the package.

To move on and conduct your own analysis, you will need access to
the LINCS Phase I Level 3 data file.  This file may already be available from your local bioinformatics core.  If you need or want to download it yourself, note that the datafile is quite large (~40GB). A robust multithreaded download client makes this much faster and less prone to failures due to connectivity hiccups.  For example, you might try:

```{r, engine = 'bash', eval=FALSE}
aria2c -x 8 -s 8 https://goo.gl/3TigFI
gzip -d *.gz
```

Alternatively, you can use the slinky library itself to fetch this file.  Again, this will take a while and may fail if your connection is not rock solid.

```{r, echo=TRUE, message=F, warning=F, eval=FALSE}
sl <- Slinky$new(key="ignore")
sl$download(type = "expression", level="3")
```

You will likely also want the metadata describing the instances in that file.  You could download it yourself, but note that the slinky package will fetch it automatically if and when it needs it, and even try to place it in the package installation directory so it is there the next time. (If it cannot be moved to that directory due to permission issues, it will just be left in your current working directory.  The `slinky` package also checks in the current working directory for the file before downloading it again.)

```{r, echo=TRUE, message=F, warning=F, eval=FALSE}
sl$download(type = "info")
```


Details on the LINCS data files are available here: https://docs.google.com/document/d/1q2gciWRhVCAAnlvF2iRLuJ7whrGP6QjpsCMq1yWz7dU

## API Key

Access to the LINCS api and clue.io requires a user key, which can be obtained from the LINCS by registering at https://clue.io (free for non-commercial use).  You can provide your key directly in the call to `Slinky$new("your_key_here")`, or define the environment variable `CLUE_API_KEY` and set it to your key.


```{r, echo=TRUE, message=F, warning=F, eval=FALSE}
library(slinky)

#update following lines with your details:
key <- "YOUR_API_KEY"
gctx <- "/path/to/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"
sl <- Slinky$new(key, gctx)

```

Alternatively, you can store your key as a single line in a text file so it does not get included in your project files (which may, among other thing, be stored in public repositories, etc.)

```{r, echo=FALSE, message=F, warning=F, eval=FALSE}
library(slinky)

#update following lines with your details:
key <- readLines("/path/to/YOUR_API_KEY")
gctx <- "/path/to/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"
sl <- Slinky$new(key, gctx)

```

# Quick Start

```{r, eval=FALSE}

library(slinky)
sl <- Slinky$new(key, 
                 system.file("extdata", "demo.gctx", package="slinky"),
                 system.file("extdata", "demo_inst_info.txt", package="slinky"))
amox_gold <- sl$clue.instances(where_clause=list("pert_type"="trt_cp",
                 "pert_iname"="amoxicillin",
                 "cell_id" = "MCF7",
                 "is_gold"=TRUE), poscon = "omit")               
amox_gold_sumexp <- sl$toSummarizedExperiment(ids = amox_gold)

```

# Further reading

See the vignette in this package for further documentation and examples.
