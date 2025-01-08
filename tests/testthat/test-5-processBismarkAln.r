context("processBismarkAln checks")

test_that("expect that there is an error when reading an unsorted sam file", {
  expect_that(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.unsorted.min.sam",
                             package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      verbose = FALSE
    ) ,
    throws_error()
  )
  expect_that(
    processBismarkAln(
      location = system.file(
        "extdata",
        "test.fastq_bismark.unsorted_chr.min.sam",
        package = "methylKit"
      ),
      sample.id = "test1",
      assembly = "hg18",
      verbose = FALSE
    ) ,
    throws_error()
  )
})



test_that("check that bismark single-end bam can be read",
          {
            expect_is(
              suppressWarnings(
              processBismarkAln(
                location = system.file("extdata", "test.bismark_single_end.sorted.bam", package = "methylKit"),
                sample.id = "test1",
                assembly = "hg18",
                treatment = NULL,
                verbose = FALSE
              )),
              'methylRaw'
            )
          })

test_that("check that bismark paired-end bam  can be read",
          {
            expect_is(suppressWarnings(
              processBismarkAln(
                location = system.file("extdata", "ctrl.bismark_paired_end.sorted.bam", package = "methylKit"),
                sample.id = "ctrl1",
                assembly = "hg18",
                treatment = NULL,
                verbose = FALSE
            )),
            'methylRaw')
          })

test_that("check that sam file can be read-in", {
  expect_is(suppressMessages(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      verbose = FALSE
    )
  ),
  'methylRaw')
  expect_is(
    suppressMessages(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      save.folder = "methylDB",
      save.db = TRUE,
      verbose = FALSE
  )),
  'methylRawDB')
})

test_that("check that reading sam and bam returns the same object", {
  expect_equal(
    processBismarkAln(
      location = system.file("extdata",
        "ctrl.bismark_paired_end.sorted.bam",
        package = "methylKit"
      ),
      sample.id = "ctrl1",
      assembly = "hg18",
      verbose = FALSE,
      mincov = 1
    ),
    processBismarkAln(
      location = system.file("extdata",
        "ctrl.bismark_paired_end.sorted.wh.sam",
        package = "methylKit"
      ),
      sample.id = "ctrl1",
      assembly = "hg18",
      verbose = FALSE,
      mincov = 1
    )
  )
  expect_equal(
    processBismarkAln(
      location = system.file("extdata",
        "test.bismark_single_end.sorted.bam",
        package = "methylKit"
      ),
      sample.id = "test",
      assembly = "hg18",
      verbose = FALSE,
      mincov = 1
    ),
    processBismarkAln(
      location = system.file("extdata",
        "test.bismark_single_end.sorted.wh.sam",
        package = "methylKit"
      ),
      sample.id = "test",
      assembly = "hg18",
      verbose = FALSE,
      mincov = 1
    )
  )
})

test_that("check that CHG context can be read in", {
  expect_is(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      read.context = "CHG",
      verbose = FALSE
  ),
  'methylRaw')
})

test_that("check that CHG context can be read-in as methylrawdb", {
  expect_is(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      read.context = "CHG",
      save.folder = "methylDB",
      save.db = TRUE,
      verbose = FALSE
  ),
  'methylRawDB')
})

test_that("check that CHH context can be read-inas methylraw", {
  expect_is(suppressMessages(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      read.context = "CHH",
      verbose = FALSE
    )
  ),
  'methylRaw')
})

test_that("check that CHH context can be read-in as methylrawdb", {
  expect_is(suppressWarnings(
    processBismarkAln(
      location = system.file("extdata", "test.fastq_bismark.sorted.min.sam", package = "methylKit"),
      sample.id = "test1",
      assembly = "hg18",
      read.context = "CHH",
      save.folder = "methylDB",
      save.db = TRUE,
      verbose = FALSE
    )
  ),
  'methylRawDB')
})

# reading multiple files
file.list2=list(system.file("extdata", "test.fastq_bismark.sorted.min.sam",
                            package = "methylKit"),
                system.file("extdata", "test.fastq_bismark.sorted.min.sam",
                            package = "methylKit"),
                system.file("extdata", "test.fastq_bismark.sorted.min.sam",
                            package = "methylKit"),
                system.file("extdata", "test.fastq_bismark.sorted.min.sam",
                            package = "methylKit"))

test_that("check that multiple files can be read-in as methylrawList", {
  expect_is(
    processBismarkAln(
      location = file.list2,
      sample.id = list("test1", "test2", "ctrl1", "ctrl1"),
      assembly = "hg18",
      save.folder = NULL,
      save.context = NULL,
      read.context = "CpG",
      nolap = FALSE,
      mincov = 10,
      minqual = 20,
      phred64 = FALSE,
      treatment = c(1, 1, 0, 0),
      verbose = FALSE
    ),
    'methylRawList'
  )
})

# test for consistent output (see https://github.com/al2na/methylKit/issues/334)
aln_sam <- 
'@HD	VN:1.0	SO:coordinate
@SQ	SN:chr19	LN:7000000
@PG	ID:Bismark	VN:v0.24.2	CL:"bismark --non_directional --parallel 8 --genome Big_Data/Test_WGBS/mapping_BOWTIE2/bismark_genome --bowtie2 --score_min L,0,-0.6 -N 0 --bam --gzip --output_dir Big_Data/Test_WGBS/mapping_BOWTIE2/bam_byName/ None -1 Big_Data/Test_WGBS/trimmed/SRR11806589_sub500000_chr19_R1_clean.fastq.gz -2 Big_Data/Test_WGBS/trimmed/SRR11806589_sub500000_chr19_R2_clean.fastq.gz"
@PG	ID:samtools	PN:samtools	PP:Bismark	VN:1.20	CL:/opt/conda/envs/methylator_v1.1_noversion/bin/samtools view -bSh -
@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.20	CL:/opt/conda/envs/methylator_v1.1_noversion/bin/samtools view -h /shared/projects/bi4edc/WGBSflow/Big_Data/Test_WGBS/mapping_BOWTIE2/bam_byName/SRR11806589_sub500000_chr19_R1_clean.fastq.gz.temp.1.gz_bismark_bt2_pe.bam
@PG	ID:samtools.2	PN:samtools	PP:samtools.1	VN:1.20	CL:/opt/conda/envs/methylator_v1.1_noversion/bin/samtools view -bSh -
@PG	ID:samtools.3	PN:samtools	PP:samtools.2	VN:1.20	CL:/opt/conda/envs/methylator_v1.1_noversion/bin/samtools view -h --threads 1 Big_Data/Test_WGBS/mapping_BOWTIE2/bam_byName/SRR11806589_sub500000_chr19_clean_bismark_bt2_pe.bam
@PG	ID:samtools.4	PN:samtools	PP:samtools.3	VN:1.20	CL:/opt/conda/envs/methylator_v1.1_noversion/bin/samtools view -bS --threads 1 -
@PG	ID:samtools.5	PN:samtools	PP:samtools.4	VN:1.20	CL:samtools sort -o Big_Data/Test_WGBS/mapping_BOWTIE2/Deduplicate/SRR11806589_sub500000_chr19.deduplicated.sort.bam -T tmp_SRR11806589_sub500000_chr19/ --threads 8 Big_Data/Test_WGBS/mapping_BOWTIE2/Deduplicate/SRR11806589_sub500000_chr19.deduplicated.bam
@PG	ID:samtools.6	PN:samtools	PP:samtools.5	VN:1.14	CL:samtools view -bh test/SRR11806589_sub500000_chr19.deduplicated.sort.bam chr19:3079455-3079616
@PG	ID:samtools.7	PN:samtools	PP:samtools.6	VN:1.14	CL:samtools view -h reduced.bam
SRR11806589.332334077_332334077_length=100	147	chr19	3079387	7	67M6I6M	=	3079292	-168	TTTAGGAGAGGTATTGTTGGATTTTTTGGTAGTATTATGTTTAGTTTTTTGAGGAATTGTTAGATTGATTTTTATCTTC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFF:FF	NM:i:25	MD:Z:0C0C0C13C5C0C1C0C7C5C0C6C7C0C1C0C3C4T1C1	XM:Z:hhx.............x.....hh.xz.......h.....hx......x.......xz.hx...x............hX	XR:Z:GA	XG:Z:CT
SRR11806589.332350504_332350504_length=100	163	chr19	3079399	6	81M	=	3079429	111	ATTACTAAATCCTCCAATAATACTATATCCAATTTTCTAAAAAACCACCAAACTAATTTCCAAAATAATTATACAAACCTA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFF	NM:i:21	MD:Z:3G2G0G7G0G2G6G4G6G1G0G4G3G3G7G1G1G0G2G5G3G0	XM:Z:...h..xh.......zx..h......h....x......x.hh....z...x...x.......x.h.hh..h.....h...x	XR:Z:GA	XG:Z:GA
SRR11806589.332355738_332355738_length=100	147	chr19	3079400	7	81M	=	3079269	-212	TTGTTGGATTTTTTGGTAGTATTATGTTTAGTTTTTTGAGGAATTGTTAGATTGATTTTTAGAGTGGTTGTATAAGTCTCC	FFFFFFFFF:F:FFFF:FF:FFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFF	NM:i:19	MD:Z:3C5C0C1C0C7C5C0C6C7C0C1C0C3C6C0C12C3C2G1	XM:Z:...x.....hh.xz.......h.....hx......x.......xz.hx...x......hx............h...hX..H	XR:Z:GA	XG:Z:CT
SRR11806589.332362696_332362696_length=100	163	chr19	3079416	6	81M	=	3079510	175	TAATACTATATCCAATTTTCTAAAAAACCACCAAACTAATTTCCAAAATAATTATACAAACCTACACTCCCACCAACAATA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFF:FF:F:FFFF,F	NM:i:17	MD:Z:2G6G4G6G1G0G4G3G3G7G1G1G0G2G5G3G16G0	XM:Z:..h......h....x......x.hh....z...x...x.......x.h.hh..h.....h...x................h	XR:Z:GA	XG:Z:GA
SRR11806589.332370546_332370546_length=100	147	chr19	3079418	6	80M	=	3079337	-161	GTATTATGTTTAGTTTTTTGAGGAATTGTTAGATTGATTTTTAGAGTGGTTGTATAAGTTTGTATTTTTATTAATAATGG	FFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:22	MD:Z:3C5C0C6C7C0C1C0C3C6C0C12C3C0C2C1C1C0C0C1C0C2C5	XM:Z:...h.....hx......x.......xz.hx...x......hx............h...hx..h.h.hhh.hh..h.....	XR:Z:GA	XG:Z:CT
SRR11806589.216596805_216596805_length=100	163	chr19	3079419	31	65M	=	3079506	153	TACTATATCCAATTTTCTAAAAAACCGCCAAACTAATTTCCAAAATAATTATACAAACCTACACT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF	NM:i:14	MD:Z:6G4G6G1G0G8G3G7G1G1G0G2G5G3G4	XM:Z:......h....x......x.hh....Z...x...x.......x.h.hh..h.....h...x....	XR:Z:GA	XG:Z:GA
SRR11806589.332350504_332350504_length=100	83	chr19	3079429	6	81M	=	3079399	-111	AATTTTCTAAAAAACCACCAAACTAATTTCCAAAATAATTATACAAACCTACACTCCCACCAACAATAAAAAAATATTCCT	FFFFFFFFFFFFFFFFFFFF:,FFF,FFFFFFF:FFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FF:FF:FFFF	NM:i:20	MD:Z:1G6G1G0G4G3G3G7G1G1G0G2G5G3G16G0G1G0G1G1G5	XM:Z:.x......x.hh....z...x...x.......x.h.hh..h.....h...x................hh.hh.h.h.....	XR:Z:CT	XG:Z:GA
SRR11806589.332373229_332373229_length=100	163	chr19	3079442	6	81M	=	3079555	192	ACCACCAAACTAATTTCCAAAATAATTATACAAACCTACACTCCCACCAACAATAAAAAAATATTCCTCTTTCTCCACATC	FFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFFFFFF	NM:i:16	MD:Z:3G3G3G7G1G1G0G2G5G3G16G0G1G0G1G1G18	XM:Z:...z...x...x.......x.h.hh..h.....h...x................hh.hh.h.h..................	XR:Z:GA	XG:Z:GA
SRR11806589.332360085_332360085_length=100	99	chr19	3079448	12	67M	=	3079569	201	AGATTGATTTTTAGAGTGGTTGTATAAGTTTGTATTTTTATTAGTAATGGAGGAGTGTTTTTTTTTT	FFFFFFFF:FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:19	MD:Z:3C6C0C12C3C0C2C1C1C0C0C1C0C1A0C14C0C1C3C0	XM:Z:...x......hx............h...hx..h.h.hhh.hh..h..............hh.h...h	XR:Z:CT	XG:Z:CT
SRR11806589.332355734_332355734_length=100	99	chr19	3079450	7	55M1D26M	=	3079577	207	ATTGATTTTTAGAGTGGTTGTATAAGTTTGTATTTTTATTAATAATGGAGGAGTGTTTTTTTTTTTTATATTTTCGTTAGT	FF:FFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFF,:F:F,F	NM:i:27	MD:Z:1C6C0C12C3C0C2C1C1C0C0C1C0C2C12^T1C0C1C3C1C0C1C2C0C3A0C2C0	XM:Z:.x......hx............h...hx..h.h.hhh.hh..h.............hh.h...h.hh.h..hh.Z..x..h	XR:Z:CT	XG:Z:CT
SRR11806589.332369609_332369609_length=100	163	chr19	3079467	7	41M2D40M	=	3079525	137	TTATACAAACCTACAATCCCACCAACAATAAAAAAATATTCCTTTCTCCACATCCTCAACAACATCTACTATCACCTAAAT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:17	MD:Z:2G5G3G2C13G0G1G0G1G1G3^CT16G3G5G2G6G3	XM:Z:..h.....h...x................hh.hh.h.h...................z...x.....x..x......x...	XR:Z:GA	XG:Z:GA
SRR11806589.332334753_332334753_length=100	99	chr19	3079483	6	22M1I58M	=	3079639	237	TTTTATTAATAATGGAGGAGTGTTTTTTTTTTTTTTATATTTTTGTTAGTATTTGTTGTTATTTGAATTTTTGATTTTAGT	FFFFFFFFFFFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFFFFFFFFFFFFF,,	NM:i:27	MD:Z:1C0C0C1C0C2C14C0C1C3C1C0C1C2C0C1C1A0C2C2C2C3C1C0C12C4C0	XM:Z:.hhh.hh..h...............hh.h...h.hh.h..hh.z..x..h..x..x...h.hx............h....h	XR:Z:CT	XG:Z:CT
SRR11806589.332331913_332331913_length=100	99	chr19	3079502	6	80M	=	3079563	142	GTGTTTTTTTTTTTTTATATTTTTGATAGTATTTGTTGTTATTTGAATTTTTGATTTTAGTTATTTTGATTGGTGTGAGG	FFFFFFFFFFFFFFFFFFF,FFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:FFFFFFFFFFF:FFFFF:FFFFF	NM:i:22	MD:Z:5C0C1C3C1C0C1C2C0C1C2C2C2C2C3C1C0C12C4C0C3C3C10	XM:Z:.....hh.h...h.hh.h..hh.z..x..h..x..x...h.hx............h....hh...x...x..........	XR:Z:CT	XG:Z:CT
SRR11806589.216596805_216596805_length=100	83	chr19	3079506	31	66M	=	3079419	-153	TCCTCTTTCTCCACATCCTCGACAACATCTACTATCACCTAAATTTTTAATCTTAACCATTCTAAC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF:::F,FFFFF:F	NM:i:7	MD:Z:24G5G2G6G7G6G7G2	XM:Z:....................Z...x.....x..x......x.......h......h.......x..	XR:Z:CT	XG:Z:GA
SRR11806589.332362696_332362696_length=100	83	chr19	3079510	6	81M	=	3079416	-175	CTTTCTCCACATCCTCACCAACATCTACTATCACCTAAATTTTTAATCTTAACCATTCTAACTAATATAAAATAAAATCTC	FFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF	NM:i:17	MD:Z:16G0A2G5G2G6G7G6G7G3G0G1G1G1G0G1G0G6	XM:Z:................z...x.....x..x......x.......h......h.......x...xh.h.h.hh.hh......	XR:Z:CT	XG:Z:GA
SRR11806589.332369609_332369609_length=100	83	chr19	3079525	7	79M	=	3079467	-137	CAACAACATCTACTATCACCTAAATTTTTAATCTTAACCATTCTAACTAATATAAAATAAAATCTCAAAATTATTTTAA	FF:FFFFFFFFFFF:FFFFFFFFFFFF,FFFF::FFFFFFFFF:F:FFFF:FFFFFFFFF:FFFFFFFFFFFFFFFFFF	NM:i:21	MD:Z:1G3G5G2G6G7G6G7G3G0G1G1G1G0G1G0G7G0G0G2G4G1	XM:Z:.z...x.....x..x......x.......h......h.......x...xh.h.h.hh.hh.......xhh..h....h.	XR:Z:CT	XG:Z:GA
SRR11806589.332335201_332335201_length=100	163	chr19	3079536	6	80M	=	3079610	154	ACTATCACCTAAATTTTTAATCTTAACCATTCTAACTAATATAAAATAAAATCTCAAAATTATTTTAATTTACATTTCCC	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFF	NM:i:20	MD:Z:0G2G6G7G6G7G3G0G1G1G1G0G1G0G7G0G0G2G4G4G8	XM:Z:x..x......x.......h......h.......x...xh.h.h.hh.hh.......xhh..h....h....h........	XR:Z:GA	XG:Z:GA
SRR11806589.332359309_332359309_length=100	163	chr19	3079551	2	81M	=	3079605	134	GAGAATCTTAACCATTCTAACTAATATAAAATAAAATCTCAAAATTATTTTAATTTACATTTCCCTATTAATTAAAAATAT	FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:26	MD:Z:0T0T0T0G6G7G3G0G1G1G1G0G1G0G7G0G0G2G4G4G9G0A1G5G0G2G1	XM:Z:...h......h.......x...xh.h.h.hh.hh.......xhh..h....h....h.........x..h.....hh..h.	XR:Z:GA	XG:Z:GA
SRR11806589.332373229_332373229_length=100	83	chr19	3079555	6	79M	=	3079442	-192	ATCTTAACCATTCTAACTAATATAAAATAAAATCTCAAAATTATTTTAATTTACATTTCCCTAATAATTAAAAATATTA	FFFFFFFFFFF:FFF,FFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF	NM:i:22	MD:Z:6G7G3G0G1G1G1G0G1G0G7G0G0G2G4G4G9G2G5G0G2G2G0	XM:Z:......h.......x...xh.h.h.hh.hh.......xhh..h....h....h.........x..h.....hh..h..h	XR:Z:CT	XG:Z:GA
SRR11806589.332374487_332374487_length=100	163	chr19	3079559	2	38M1D42M	=	3079653	174	TCCCCATTCTAACTAATATAAAATAAAATCTCAATATTTTTTAATTTACATTTCCCTAATAATTAAAAATACTAAACATT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF:F	NM:i:24	MD:Z:1A0G7G3G0G1G1G1G0G1G0G7G0G0G2^G4G4G9G2G5G0G2G0T1G6	XM:Z:..........x...xh.h.h.hh.hh.......x.h......h....h.........x..h.....hh..h..h......	XR:Z:GA	XG:Z:GA
SRR11806589.332331913_332331913_length=100	147	chr19	3079563	6	81M	=	3079502	-142	TATTTTGATTGGTGTGAGGTGGAATTTTAGGGTTGTTTTGATTTGTATTTTATTGATGATTAAGGATGTTGAATATTTTTT	FFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFF	NM:i:10	MD:Z:0C3C3C16C1C17C4C0C0C20C7	XM:Z:h...x...x................h.x.................h....h.x....................h.......	XR:Z:GA	XG:Z:CT
SRR11806589.332360085_332360085_length=100	147	chr19	3079569	12	80M	=	3079448	-201	GATTGGTGTGAGGTGGAATTTTAGGGTTGTTTTGATTTGTATTTTTTTGATGATTAAGGATGTTGAATATTTTTTTAAGT	FFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:9	MD:Z:2C16C1C17C4C0C0C20C7C4	XM:Z:..x................h.x.................h....hhx....................h.......h....	XR:Z:GA	XG:Z:CT
SRR11806589.332355734_332355734_length=100	147	chr19	3079577	7	80M	=	3079450	-207	TGAGGTGGAATATTAGGGTTGTTTTGATTTGTATTTTTTTGATGATTAAGGATGTTGAATATTTTTTTAAGTGTTTTTTT	FFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:11	MD:Z:11C1C17C4C0C0C20C7C5C2C1C1	XM:Z:.............x.................h....hhx....................h.......h.....h..h.x.	XR:Z:GA	XG:Z:CT
SRR11806589.332348899_332348899_length=100	99	chr19	3079600	6	80M	=	3079656	135	TTGATTTGTATTTTTTTGATGATTAAGGATTTTGAGTATTTTTTTAAGTGTTTTTTTGTTATTTGGTATTTTTTAGGTGA	FF,FFFFFFFF,:FFFFFFF,FFFFF:FFF,FF,FFFFFFF,F:FFFF:FFFFFFFF:FFF::FFFFFFFF:FFFF,:FF	NM:i:17	MD:Z:8C4C0C0C14G4A0C7C5C2C1C2C0C3C6C0C1C6	XM:Z:........h....hhx....................h.......h.....h..h.x..hh...z......hh.x......	XR:Z:CT	XG:Z:CT
SRR11806589.332359309_332359309_length=100	83	chr19	3079605	2	80M	=	3079551	-134	TTACATTTCCCTATTAATTAAAAATATTAAACATTTTTTCAAATACTTCTCTACCATTCGATATTACTCAAATAAAAATT	FFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	NM:i:17	MD:Z:2G9G0A1G5G0G2G2G13G1G7G7G4C4G0G1G1G4	XM:Z:..h.........x..h.....hh..h..h.............h.h.......x......Zx.........xh.h.h....	XR:Z:CT	XG:Z:GA
SRR11806589.332356300_332356300_length=100	163	chr19	3079610	2	80M	=	3079795	261	TTTCCCTAATGATTAAAAATACTAAACATTTTTTCAAATACTTCTCTACCATTCGATATTCCTCAAATAAAAATTCTTCG	:F:FF:FFF:FFF,FFF,F:FFFF:FFF,FF,FFFFFFFF,FF:F:F:,FF:FFFFFF,FFF:FFFF:FFFFF,::F:FF	NM:i:15	MD:Z:7G8G0G2G0T1G13G1G7G7G9G0G1G1G7T1	XM:Z:.......x..H.....hh..h..h.............h.h.......x......Zx.........xh.h.h........H	XR:Z:GA	XG:Z:GA
SRR11806589.332335201_332335201_length=100	83	chr19	3079610	6	80M	=	3079536	-154	TTTCCCTAATAATTAAAAATATTAAACATTTTTTCAAATACTTCTCTACCATTCAATATTCCTCAAATAAAAATTCTTTA	FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFF	NM:i:16	MD:Z:7G2G5G0G2G2G13G1G7G6G0G9G0G1G1G8G0	XM:Z:.......x..h.....hh..h..h.............h.h.......x......zx.........xh.h.h........h	XR:Z:CT	XG:Z:GA'

writeLines(aln_sam,con = "aln.sam")

aln_out <- 
'"chr","start","end","strand","coverage","numCs","numTs"
"chr19",3079413,3079413,"+",2,0,2
"chr19",3079414,3079414,"-",1,0,1
"chr19",3079444,3079444,"+",3,0,3
"chr19",3079445,3079445,"-",4,1,3
"chr19",3079525,3079525,"+",2,0,2
"chr19",3079526,3079526,"-",4,1,3
"chr19",3079663,3079663,"+",1,0,1
"chr19",3079664,3079664,"-",3,2,1'

aln_df <- read.csv(text = aln_out,header = TRUE, sep = ",")

test_that("check that output is consisten with methylKit.1.20.0", {
  expect_identical(
    getData(
      processBismarkAln(
        location ="aln.sam",
        sample.id = "test",
        assembly = "mm39",
        save.folder = NULL,
        save.context = NULL,
        read.context = "CpG",
        nolap = FALSE,
        mincov = 1,
        minqual = 20,
        phred64 = FALSE,
        treatment = 0,
        verbose = FALSE
      )
    ),
    aln_df
  )
})

unlink("aln.sam")
unlink("tests/testthat/methylDB",recursive = TRUE)