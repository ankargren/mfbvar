.onUnload <- function (libpath) {
  library.dynam.unload("mfbvar", libpath)
}
