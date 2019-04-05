context("Tidy Up")

dirs <- list.dirs(path = ".",full.names = TRUE,recursive = TRUE)
dirs <- grep(pattern = "methylDB",x = dirs,value = TRUE)
unlink(x = dirs, recursive = TRUE)
unlink("Rplots.pdf")
unlink("test.bed")
