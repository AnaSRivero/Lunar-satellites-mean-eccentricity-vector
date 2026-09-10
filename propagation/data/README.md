# Data required

This code needs two things placed in one local folder (referred to as
`orekitDataDir` in the scripts):

1. **The standard Orekit physical data set**, used for time scales, frames,
   ephemerides, etc.:
   https://gitlab.orekit.org/orekit/orekit-data
   (clone/download and unzip as described in that repository's own README).

2. **The lunar gravity field, GRGM1200B**, provided here as `lunar_grail.gfc`.
   This is a public NASA GSFC/JPL model (based on GRAIL data), originally
   distributed via the ICGEM service:
   http://icgem.gfz-potsdam.de/tom_longtime

   The model is published up to degree/order 1200; this study only uses
   truncations up to degree/order 50 (see `getMoonPropagator.m`). Copy
   this file into the same folder as the `orekit-data` set above.

   Note: the header of the original ICGEM file was missing the
   `tide_system`, `errors` and `norm` keys, which some Orekit versions
   (e.g. 13.1.4) require and otherwise reject with an "unexpected file
   format" error. The header shipped here adds those keys (values
   `tide_free` / `formal` / `fully_normalized`, matching the model's
   actual properties); the coefficient data itself is untouched.

No other data files (space weather, DE44x ephemerides, etc.) beyond what
`orekit-data` already provides are used by this propagation code.
