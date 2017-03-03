.onUnload <- function (libpath) {
  library.dynam.unload("methylKit", libpath)
}