# Abstract

The L1000 data set from the NIH LINCS program holds the promise to deconvolute a wide range of biological questions in transcriptional 
space.  However, using this large and decentralized data set presents its own challenge.  Here we demonstrate how to use the slinky 
package to simplify filtering and accessing this data in an efficient analysis pipeline.

# Prerequisites

## Memory

Slinky uses the rhdf5 package to access slices of the L1000 data from disk.  As a result, you only need enough memory to hold the data slices you are working with, and possibly the metadata as well.  (The entire metadata data frame for all 1.3 million instances requires a modest 220MB of RAM). The development machine on which the slinky package was primarily developed has 16GB of RAM and can load and manipulate tens of thousands of instances from the L1000 dataset without difficulty.  Of course, YMMV depending on the specifications of your computer.

## Data

In this vignette we perform some basic analyses using the L1000 data expression data from the LINCS project.  To complete this vignette you will need access to the LINCS Phase I Level 3 data file.  This file may already be available from your local bioinformatics core.  If you need or want to download it yourself, note that the datafile is quite large (~40GB). A robust multithreaded download client makes this much faster and less prone to failures due to connectivity hiccups.  For example, you might try:

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


# Data selection

There are several ways to identify a useful slice of the L1000 data set to load.  The first way is to figure out the column and row indices you are interested in and specify these directly.  For example, you might do this by filtering the info file that accompanies the L1000 data from GEO.  This approach does not even require an API key since it does not call the API.

(It may be usefuly to note that the "landmark" genes are the first 978 rows of the dataset.  For our demonstrations here we will stick to that subset of genes.  The remaining rows contain the inferred expression for the rest of the transcriptome.)

## From the info file

The `loadInfo` method will find the info file (and download it if it cannot be located) and then load it into the `metadata` slot of your Slinky object.

```{r, echo=FALSE, message=F, warning=F, eval=FALSE}
gctx <- "/path/to/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"
sl <- Slinky$new("none", gctx)
sl$loadInfo()

col.ix <- which(sl$metadata$pert_iname == "amoxicillin" & sl$metadata$cell_id == "MCF7")
data <- sl$readGCTX(index = list(1:978, col.ix))
```

Note that the `index` argument to the `readGCTX` expects a list (of row indices followed by column indices), because that is what `rhdf5::h5read` requires.

Note that loading contiguous data out of the gctx file is much faster than data spread throughout the file.  For example, loading 1000 consecutive instances (over the network) takes 3.4 seconds on my machine, while loading 1000 random instances takes almost 10 times that long.  Unfortunately, the order of the instances in the gctx file is different than in the info file.  However, you can retrieve the complete list of instances from the gctx file (in order) using the `colnames` method (e.g. `ids <- sl$colnames()`). Depending on your use case, you may be able to leverage this information to optimize performance.


## By querying the API

While it is not as fast, querying the clue.io API is simpler, and also allows you to query on metadata not present in the info file provided via GEO.  For example, we could identify all the instances that are not only perturbed by amoxicillin, but are also designated as "gold" by the LINCS program (by virtue of having low variance across replicates):

```{r, echo=TRUE, message=F, warning=F, eval=FALSE}
amox_gold <- sl$clue.instances(where_clause=list("pert_type"="trt_cp",
                                                "pert_iname"="amoxicillin",
                                                "is_gold"=TRUE), poscon = "omit")
```


Note the `poscon = "omit"` argument.  There are a small number of instances that have a `pert_type` of `trt_poscon` in clue.io's `profiles` endpoint but are recoded as `trt_cp` in the `sigs` endpoint.  This can cause confusion/errors in downstream analysis, so by these samples can be omitted from query results.  In fact, this is the default behavior of the `clue.instances` method (the argument was simply specified above for emphasis).

You could then load these instances by "brute force".

```{r, echo=TRUE, message=F, warning=F, eval=FALSE}
ix <- which sl$colnames() %in% amox_gold
amox_gold_data <- sl$readGCTX(index = list(1:978, ix))
```

But it is usually more useful to keep this data annotated with its corresponding metadata in the form of an `ExpressionSet`.  This can be easily done with the `toEset` method.


```{r, echo=TRUE, message=F, warning=F, eval=FALSE}

amox_gold_eset <- sl$toEset(index = list( 1:978, ix))

```

Or you can specify the ids directly.

```{r, echo=TRUE, message=F, warning=F, eval=FALSE}

amox_gold_eset <- sl$toEset(ids = amox_gold)

```

In fact, you could have skipped the trouble of looking up the ids in the first place and simply passed the `where_clause` to `toEset` directly.

```{r, echo=TRUE, message=F, warning=F, eval=FALSE}

amox_gold_eset <- sl$toEset(where_clause=list("pert_type"="trt_cp",
                                                "pert_iname"="amoxicillin",
                                                "is_gold"=TRUE))

```

Note that by default, toEset only loads the L1000 genes, so there is no need to specify the row index.  If you desire all the genes, simply specify `inferred=TRUE` when calling `toEset`.  If you want some other subset of genes, you will need to use the `index` method shown previously, and not the ids method. This underscores the fact that `toEset` expects exactly one of the `ids`, `index`, or `where_clause` argument.  Specifying more than one of those arguments, or none of them, will throw an error.


Here is a more complex example.  Let us identify those instances treated with FDA approved compounds that are part of the "gold" subset of highly consistent instances (as defined by LINCS).  Definitively identifying FDA approved drugs is a little tricky due to subtle (and not so subtle) differences in vocabularies, but we exploit the "repositioning" data available at the `rep_drugs` endpoint to achieve our ends. 


```{r, echo=FALSE, message=F, warning=F, eval=FALSE}


fda <- sl$clue("rep_drugs", where_clause=list("status_source"=list(like = "FDA Orange"),
                                                   "final_status"="Launched",
                                                   "animal_only"="0",
                                                   "in_cmap"=TRUE),
                                                   verbose=TRUE)

fda_pert <- sl$clue.instances(where_clause=list("pert_type"="trt_cp", 
                                                "is_gold"=TRUE,
                                                "pert_iname"=list("inq"=fda$pert_iname)),
                              verbose=TRUE)

```

Note the `inq` syntax used above to provide an array of possible values to match.  This is a vernacular of the `loopback` framework on which the clue.io service is built.

# Further reading

See the vignette in this package for further documntation and examples.
