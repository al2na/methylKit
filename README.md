methylKit
========

This is where methylKit project will be developed from 2014. The [google code homepage](https://code.google.com/p/methylkit/) 
will be active for the foreseeable future.


Installation
---------
```R 
# dependencies
install.packages( c("data.table","devtools"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges","IRanges"))

# install from github
library(devtools)
install_github("methylKit", username = "al2na")
```

Citing methylKit
---------

If you used methylKit please cite:


 * Altuna Akalin, Matthias Kormaksson, Sheng Li, Francine E. Garrett-Bakelman, Maria E. Figueroa, Ari Melnick, Christopher E. Mason. _(2012)_. *"[http://genomebiology.com/2012/13/10/R87/ methylKit: A comprehensive R package for the analysis of genome-wide DNA methylation profiles.]"* _Genome Biology_ , 13:R87.


and also consider citing the following publication as a use-case with specific cutoffs:

 * Altuna Akalin, Francine E. Garrett-Bakelman, Matthias Kormaksson, Jennifer Busuttil, Lu Zhang, Irina Khrebtukova, Thomas A. Milne, Yongsheng Huang, Debabrata Biswas, Jay L. Hess, C. David Allis, Robert G. Roeder, Peter J. M. Valk, Bob Löwenberg, Ruud Delwel, Hugo F. Fernandez, Elisabeth Paietta, Martin S. Tallman, Gary P. Schroth, Christopher E. Mason, Ari Melnick, Maria E. Figueroa. _(2012)_. *"[http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002781 Base-Pair Resolution DNA Methylation Sequencing Reveals Profoundly Divergent Epigenetic Landscapes in Acute Myeloid Leukemia.]"* _PLoS Genetics_ 8(6).


Contact & Questions
-------
e-mail to [mailto:methylkit_discussion@googlegroups.com methylkit_discussion@googlegroups.com] or post a question using [https://groups.google.com/forum/#!forum/methylkit_discussion the web interface]

if you are going to submit bug reports or ask questions, please send sessionInfo() output from R console as well.

Questions are very welcome, although we suggest you read the paper, documentation(function help pages and [http://methylkit.googlecode.com/git/inst/doc/methylKit.pdf vignette]) and [http://zvfak.blogspot.com/search/label/methylKit blog entries] first. The answer to your question might be there already.


License
---------
Artistic License/GPL
