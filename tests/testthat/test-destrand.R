context("test for known destranding issues")

## reproduce issue 296
write(file = "ctrl_cytosineReport.tst",x = 
        "chr12__pp	3000261	-	0	11	CHH	CCT
chr12__pp	3000262	-	0	12	CHH	CCC
chr12__pp	3000264	+	0	17	CHH	CCA
chr12__pp	3000265	+	0	17	CHG	CAG
chr12__pp	3000267	-	0	13	CHG	CTG
chr12__pp	3000270	+	0	17	CHG	CCG
chr12__pp	3000271	+	17	0	CG	CGG
chr12__pp	3000286	-	0	14	CHH	CCC
chr12__pp	3000289	-	0	13	CHH	CTT
chr1__pp	3000286	-	0	14	CHH	CCC
chr1__pp	3000289	-	0	13	CHH	CTT
chr1__pp	3000290	-	0	13	CHH	CCT
chr1__pp	3000291	-	0	13	CHH	CCC
chr1__pp	3000292	-	0	12	CHH	CCC
chr1__pp	3000293	+	0	16	CHH	CTC
chr1__pp	3000295	+	15	1	CG	CGG
chr1__pp	3000296	-	11	1	CG	CGA
chr1__pp	3000297	-	0	12	CHG	CCG
chr1__pp	3000301	-	0	11	CHH	CAA
chr1__pp	3000303	-	0	11	CHH	CAC
chr1__pp	3000310	+	0	16	CHH	CAC
chr12__pp	3000290	-	0	13	CHH	CCT
chr12__pp	3000291	-	0	13	CHH	CCC
chr12__pp	3000292	-	0	12	CHH	CCC
chr12__pp	3000293	+	0	16	CHH	CTC
chr12__pp	3000295	+	15	1	CG	CGG
chr12__pp	3000296	-	11	1	CG	CGA
chr12__pp	3000297	-	0	12	CHG	CCG
chr12__pp	3000301	-	0	11	CHH	CAA
chr12__pp	3000303	-	0	11	CHH	CAC
chr12__pp	3000310	+	0	16	CHH	CAC
chr12__pp	3000312	+	0	16	CHH	CAA
chr12__pp	3000315	-	0	12	CHH	CTT
chr12__pp	3000316	-	0	12	CHH	CCT
chr1__pp0	3000261	-	0	11	CHH	CCT
chr1__pp0	3000262	-	0	12	CHH	CCC
chr1__pp0	3000264	+	0	17	CHH	CCA
chr1__pp0	3000284	-	0	14	CHH	CCA
chr1__pp0	3000285	-	0	14	CHH	CCC
chr1__pp0	3000286	-	0	14	CHH	CCC
chr1__pp0	3000289	-	0	13	CHH	CTT
chr1__pp0	3000290	-	0	13	CHH	CCT
chr1__pp0	3000291	-	0	13	CHH	CCC
chr1__pp0	3000292	-	0	12	CHH	CCC
chr1__pp0	3000293	+	0	16	CHH	CTC
chr1__pp0	3000295	+	15	1	CG	CGG
chr1__pp0	3000296	-	11	1	CG	CGA
chr1__pp0	3000297	-	0	12	CHG	CCG
chr1__pp0	3000301	-	0	11	CHH	CAA
chr1__pp0	3000303	-	0	11	CHH	CAC
chr1__pp0	3000310	+	0	16	CHH	CAC
chr1__pp0	3000312	+	0	16	CHH	CAA
chr1__pp0	3000315	-	0	12	CHH	CTT
chr1__pp0	3000316	-	0	12	CHH	CCT
chr9__pp	3000261	-	0	11	CHH	CCT
chr9__pp	3000262	-	0	12	CHH	CCC
chr9__pp	3000264	+	0	17	CHH	CCA
chr9__pp	3000265	+	0	17	CHG	CAG
chr9__pp	3000267	-	0	13	CHG	CTG
chr9__pp	3000270	+	0	17	CHG	CCG
chr9__pp	3000271	+	17	0	CG	CGG
chr9__pp	3000272	-	12	1	CG	CGG
chr9__pp	3000273	-	0	13	CHG	CCG
chr9__pp	3000276	+	0	17	CHH	CTA
chr9__pp	3000281	+	1	16	CHG	CTG
chr9__pp	3000283	-	0	13	CHG	CAG
chr9__pp	3000284	-	0	14	CHH	CCA
chr9__pp	3000285	-	0	14	CHH	CCC
chr9__pp	3000286	-	0	14	CHH	CCC
chr9__pp	3000289	-	0	13	CHH	CTT
chr9__pp	3000290	-	0	13	CHH	CCT
chr9__pp	3000291	-	0	13	CHH	CCC
chr9__pp	3000292	-	0	12	CHH	CCC
chr1__pp0	3000265	+	0	17	CHG	CAG
chr1__pp0	3000267	-	0	13	CHG	CTG
chr1__pp0	3000270	+	0	17	CHG	CCG
chr1__pp0	3000271	+	17	0	CG	CGG
chr1__pp0	3000272	-	12	1	CG	CGG
chr1__pp0	3000273	-	0	13	CHG	CCG
chr1__pp0	3000276	+	0	17	CHH	CTA
chr1__pp0	3000281	+	1	16	CHG	CTG
chr1__pp0	3000283	-	0	13	CHG	CAG
chr1__pp	3000262	-	0	12	CHH	CCC
chr1__pp	3000264	+	0	17	CHH	CCA
chr1__pp	3000265	+	0	17	CHG	CAG
chr1__pp	3000267	-	0	13	CHG	CTG
chr1__pp	3000270	+	0	17	CHG	CCG
chr1__pp	3000271	+	17	0	CG	CGG
chr1__pp	3000272	-	12	1	CG	CGG
chr1__pp	3000273	-	0	13	CHG	CCG
chr1__pp	3000276	+	0	17	CHH	CTA
chr1__pp	3000281	+	1	16	CHG	CTG
chr1__pp	3000283	-	0	13	CHG	CAG
chr1__pp	3000284	-	0	14	CHH	CCA
chr1__pp	3000285	-	0	14	CHH	CCC
chr1__pp	3000312	+	0	16	CHH	CAA
chr1__pp	3000315	-	0	12	CHH	CTT
chr1__pp	3000316	-	0	12	CHH	CCT
chr12__pp	3000272	-	12	1	CG	CGG
chr12__pp	3000273	-	0	13	CHG	CCG
chr12__pp	3000276	+	0	17	CHH	CTA
chr12__pp	3000281	+	1	16	CHG	CTG
chr12__pp	3000283	-	0	13	CHG	CAG
chr12__pp	3000284	-	0	14	CHH	CCA
chr12__pp	3000285	-	0	14	CHH	CCC
chr9__pp	3000293	+	0	16	CHH	CTC
chr9__pp	3000295	+	15	1	CG	CGG
chr9__pp	3000296	-	11	1	CG	CGA
chr9__pp	3000297	-	0	12	CHG	CCG
chr9__pp	3000301	-	0	11	CHH	CAA
chr9__pp	3000303	-	0	11	CHH	CAC
chr9__pp	3000310	+	0	16	CHH	CAC
chr9__pp	3000312	+	0	16	CHH	CAA
chr9__pp	3000315	-	0	12	CHH	CTT
chr9__pp	3000316	-	0	12	CHH	CCT
chr1__pp	3000261	-	0	11	CHH	CCT"
)

write(file="test_cytosineReport.tst", x =
        "chr1__pp2__pp	3000261	-	0	11	CHH	CCT
chr1__pp2__pp	3000262	-	0	12	CHH	CCC
chr1__pp2__pp	3000264	+	0	17	CHH	CCA
chr1__pp2__pp	3000265	+	0	17	CHG	CAG
chr1__pp2__pp	3000267	-	0	13	CHG	CTG
chr1__pp2__pp	3000270	+	0	17	CHG	CCG
chr1__pp2__pp	3000271	+	17	0	CG	CGG
chr1__pp2__pp	3000272	-	12	1	CG	CGG
chr1__pp2__pp	3000273	-	0	13	CHG	CCG
chr1__pp2__pp	3000276	+	0	17	CHH	CTA
chr1__pp2__pp	3000281	+	1	16	CHG	CTG
chr1__pp2__pp	3000283	-	0	13	CHG	CAG
chr1__pp2__pp	3000284	-	0	14	CHH	CCA
chr1__pp2__pp	3000285	-	0	14	CHH	CCC
chr1__pp2__pp	3000286	-	0	14	CHH	CCC
chr1__pp2__pp	3000289	-	0	13	CHH	CTT
chr1__pp2__pp	3000290	-	0	13	CHH	CCT
chr1__pp2__pp	3000291	-	0	13	CHH	CCC
chr1__pp2__pp	3000292	-	0	12	CHH	CCC
chr1__pp2__pp	3000293	+	0	16	CHH	CTC
chr1__pp2__pp	3000295	+	15	1	CG	CGG
chr1__pp2__pp	3000296	-	11	1	CG	CGA
chr1__pp2__pp	3000297	-	0	12	CHG	CCG
chr1__pp2__pp	3000301	-	0	11	CHH	CAA
chr1__pp2__pp	3000303	-	0	11	CHH	CAC
chr1__pp2__pp	3000310	+	0	16	CHH	CAC
chr1__pp2__pp	3000312	+	0	16	CHH	CAA
chr1__pp2__pp	3000315	-	0	12	CHH	CTT
chr1__pp2__pp	3000316	-	0	12	CHH	CCT
chr1__pp0	3000261	-	0	11	CHH	CCT
chr1__pp0	3000262	-	0	12	CHH	CCC
chr1__pp0	3000264	+	0	17	CHH	CCA
chr1__pp0	3000265	+	0	17	CHG	CAG
chr1__pp0	3000267	-	0	13	CHG	CTG
chr1__pp0	3000270	+	0	17	CHG	CCG
chr1__pp0	3000271	+	17	0	CG	CGG
chr1__pp0	3000272	-	12	1	CG	CGG
chr1__pp0	3000273	-	0	13	CHG	CCG
chr1__pp0	3000276	+	0	17	CHH	CTA
chr1__pp0	3000281	+	1	16	CHG	CTG
chr1__pp0	3000283	-	0	13	CHG	CAG
chr1__pp0	3000284	-	0	14	CHH	CCA
chr1__pp0	3000285	-	0	14	CHH	CCC
chr1__pp0	3000286	-	0	14	CHH	CCC
chr1__pp0	3000289	-	0	13	CHH	CTT
chr1__pp0	3000290	-	0	13	CHH	CCT
chr1__pp0	3000291	-	0	13	CHH	CCC
chr1__pp0	3000292	-	0	12	CHH	CCC
chr1__pp0	3000293	+	0	16	CHH	CTC
chr1__pp0	3000295	+	15	1	CG	CGG
chr1__pp0	3000296	-	11	1	CG	CGA
chr1__pp0	3000297	-	0	12	CHG	CCG
chr1__pp0	3000301	-	0	11	CHH	CAA
chr1__pp0	3000303	-	0	11	CHH	CAC
chr1__pp0	3000310	+	0	16	CHH	CAC
chr1__pp0	3000312	+	0	16	CHH	CAA
chr1__pp0	3000315	-	0	12	CHH	CTT
chr1__pp0	3000316	-	0	12	CHH	CCT
chr9__pp	3000261	-	0	11	CHH	CCT
chr9__pp	3000262	-	0	12	CHH	CCC
chr9__pp	3000264	+	0	17	CHH	CCA
chr9__pp	3000265	+	0	17	CHG	CAG
chr9__pp	3000267	-	0	13	CHG	CTG
chr9__pp	3000270	+	0	17	CHG	CCG
chr9__pp	3000271	+	17	0	CG	CGG
chr9__pp	3000272	-	12	1	CG	CGG
chr9__pp	3000273	-	0	13	CHG	CCG
chr9__pp	3000276	+	0	17	CHH	CTA
chr9__pp	3000281	+	1	16	CHG	CTG
chr9__pp	3000283	-	0	13	CHG	CAG
chr9__pp	3000284	-	0	14	CHH	CCA
chr9__pp	3000285	-	0	14	CHH	CCC
chr9__pp	3000286	-	0	14	CHH	CCC
chr9__pp	3000289	-	0	13	CHH	CTT
chr9__pp	3000290	-	0	13	CHH	CCT
chr9__pp	3000291	-	0	13	CHH	CCC
chr9__pp	3000292	-	0	12	CHH	CCC
chr9__pp	3000293	+	0	16	CHH	CTC
chr9__pp	3000295	+	15	1	CG	CGG
chr9__pp	3000296	-	11	1	CG	CGA
chr9__pp	3000297	-	0	12	CHG	CCG
chr9__pp	3000301	-	0	11	CHH	CAA
chr9__pp	3000303	-	0	11	CHH	CAC
chr9__pp	3000310	+	0	16	CHH	CAC
chr9__pp	3000312	+	0	16	CHH	CAA
chr9__pp	3000315	-	0	12	CHH	CTT
chr9__pp	3000316	-	0	12	CHH	CCT
chr1__pp	3000261	-	0	11	CHH	CCT
chr1__pp	3000262	-	0	12	CHH	CCC
chr1__pp	3000264	+	0	17	CHH	CCA
chr1__pp	3000265	+	0	17	CHG	CAG
chr1__pp	3000267	-	0	13	CHG	CTG
chr1__pp	3000270	+	0	17	CHG	CCG
chr1__pp	3000271	+	17	0	CG	CGG
chr1__pp	3000272	-	12	1	CG	CGG
chr1__pp	3000273	-	0	13	CHG	CCG
chr1__pp	3000276	+	0	17	CHH	CTA
chr1__pp	3000281	+	1	16	CHG	CTG
chr1__pp	3000283	-	0	13	CHG	CAG
chr1__pp	3000284	-	0	14	CHH	CCA
chr1__pp	3000285	-	0	14	CHH	CCC
chr1__pp	3000286	-	0	14	CHH	CCC
chr1__pp	3000289	-	0	13	CHH	CTT
chr1__pp	3000290	-	0	13	CHH	CCT
chr1__pp	3000291	-	0	13	CHH	CCC
chr1__pp	3000292	-	0	12	CHH	CCC
chr1__pp	3000293	+	0	16	CHH	CTC
chr1__pp	3000295	+	15	1	CG	CGG
chr1__pp	3000296	-	11	1	CG	CGA
chr1__pp	3000297	-	0	12	CHG	CCG
chr1__pp	3000301	-	0	11	CHH	CAA
chr1__pp	3000303	-	0	11	CHH	CAC
chr1__pp	3000310	+	0	16	CHH	CAC
chr1__pp	3000312	+	0	16	CHH	CAA
chr1__pp	3000315	-	0	12	CHH	CTT
chr1__pp	3000316	-	0	12	CHH	CCT")



files <- list("test_cytosineReport.tst","ctrl_cytosineReport.tst")

myobj_test <- methRead(files,
                       sample.id = list("test","ctrl"),
                       assembly = "misc",
                       dbtype="tabix",
                       pipeline = "bismarkCytosineReport",
                       header=TRUE,
                       sep="\t",
                       context="CpG",
                       resolution="base",
                       dbdir="CpG_methylDB_test",
                       treatment=c(2,1),
                       mincov=1)

test_that("test if unite with destranding works as tabix", {
  expect_is(meth_test_db <- unite(myobj_test, destrand=TRUE, min.per.group=1L, chunk.size = 1e1, save.db = TRUE), "methylBaseDB")
})
test_that("test if unite with destranding works in memory and as tabix", {
  expect_is(meth_test_mem <- unite(myobj_test, destrand=TRUE, min.per.group=1L, chunk.size = 1e1, save.db = FALSE), "methylBase")
})
test_that("test if output of unite with destranding is the same in memory and as tabix", {
  expect_identical(getData(meth_test_mem), getData(meth_test_db))
})
