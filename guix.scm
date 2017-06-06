;;; methylKit -- DNA methylation analysis
;;; Copyright © 2016 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of methylKit.
;;;
;;; methylKit is free software under the Artistic License 2.0.

;;; Run the following command to enter a development environment for
;;; methylKit:
;;;
;;;  $ guix environment -l guix.scm

(use-modules ((guix licenses) #:prefix license:)
             (guix packages)
             (guix download)
             (guix utils)
             (guix build-system r)
             (gnu packages)
             (gnu packages bioinformatics)
             (gnu packages statistics)
             (gnu packages perl))

(define-public r-methylkit
  (package
    (name "r-methylkit")
    (version "0.9.5")
    (source #f)
    (build-system r-build-system)
    (native-inputs
     `(("r-testthat" ,r-testthat)))
    (inputs
     `(("perl" ,perl)))
    (propagated-inputs
     `(("r-data-table" ,r-data-table)
       ("r-genomeinfodb" ,r-genomeinfodb)
       ("r-s4vectors" ,r-s4vectors)
       ("r-iranges" ,r-iranges)
       ("r-genomicranges" ,r-genomicranges)))
    (synopsis "DNA methylation analysis")
    (description
     "methylKit is an R package for DNA methylation analysis and
annotation from high-throughput bisulfite sequencing.  The package is
designed to deal with sequencing data from RRBS and its variants, but
also target-capture methods and whole genome bisulfite sequencing.  It
also has functions to analyze base-pair resolution 5hmC data from
experimental protocols such as oxBS-Seq and TAB-Seq.  Perl is needed
to read SAM files only.")
    (home-page "https://github.com/al2na/methylKit")
    (license license:artistic2.0)))

r-methylkit
