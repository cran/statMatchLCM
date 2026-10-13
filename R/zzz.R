.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "statMatchLCM: Statistical matching using latent class models.\n",
    "Start with ?datAB_to_SM or ?datABC_to_SM."
  )
}
