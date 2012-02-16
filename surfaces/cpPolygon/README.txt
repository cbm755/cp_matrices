Code by Tom Maerz, computes CP and signed distance for general polygons

This code also does tangents and various other things, these have been
(temporarily?) to "other" to keep them out of the namespace (this
directory will likely be added to matlab's path).

Issues: this code has some bugs.

TODO: BUG: signed distance depends on the orientation of the curve.

TODO: BUG: signed distance sometimes zero incorrectly (see unit tests).
