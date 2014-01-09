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


 * Altuna Akalin, Matthias Kormaksson, Sheng Li, Francine E. Garrett-Bakelman, Maria E. Figueroa, Ari Melnick, Christopher E. Mason. _(2012)_. [methylKit: A comprehensive R package for the analysis of genome-wide DNA methylation profiles.](http://genomebiology.com/2012/13/10/R87/) _Genome Biology_ , 13:R87.


and also consider citing the following publication as a use-case with specific cutoffs:

 * Altuna Akalin, Francine E. Garrett-Bakelman, Matthias Kormaksson, Jennifer Busuttil, Lu Zhang, Irina Khrebtukova, Thomas A. Milne, Yongsheng Huang, Debabrata Biswas, Jay L. Hess, C. David Allis, Robert G. Roeder, Peter J. M. Valk, Bob Löwenberg, Ruud Delwel, Hugo F. Fernandez, Elisabeth Paietta, Martin S. Tallman, Gary P. Schroth, Christopher E. Mason, Ari Melnick, Maria E. Figueroa. _(2012)_. [Base-Pair Resolution DNA Methylation Sequencing Reveals Profoundly Divergent Epigenetic Landscapes in Acute Myeloid Leukemia.](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002781) _PLoS Genetics_ 8(6).


Contact & Questions
-------
e-mail to [methylkit_discussion@googlegroups.com](mailto:methylkit_discussion@googlegroups.com) or post a question using [the web interface](https://groups.google.com/forum/#!forum/methylkit_discussion)

if you are going to submit bug reports or ask questions, please send sessionInfo() output from R console as well.

Questions are very welcome, although we suggest you read the paper, documentation(function help pages and [vignette](https://github.com/al2na/methylKit/blob/master/inst/doc/methylKit.pdf?raw=true)) and [blog entries](http://zvfak.blogspot.com/search/label/methylKit) first. The answer to your question might be there already.


Contribute to the development
---------
You can contribute to the methylKit development via github by checking out "development" branch, making your changes and doing a pull request (all of these should be done on the "development" branch NOT on the "master" branch). 
In addition:
 * Bump up the version in the DESCRIPTION file on the 4th number. For example, the master branch has the version numbering as in X.Y.Z, if you make a change to the development branch you should bump up the version to X.Y.Z.1 if there are already changes done in the development just bump up the fourth number. 
 * Add your changes to the NEWS file as well under the correct version and attribute the changes to yourself, such as "Contributed by X"

License
---------
Artistic License/GPL
