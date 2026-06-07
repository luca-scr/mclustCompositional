.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- paste("Loaded package 'mclustCompositional' version", packageVersion("mclustCompositional"))
  packageStartupMessage(msg)      
  invisible()
}
