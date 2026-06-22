# Platform sensitive path determination
psp <- function(unix, windows) if(.Platform[["OS.type"]] == "unix") unix else windows
