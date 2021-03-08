# BSpline Minimum Distance to a Point

One day, while reading literature related to B-Spline geometry, I stumbled upon [this paper] (https://hal.inria.fr/inria-00518359)

The algorithm described by the authors is both elegant and performant, and they demonstrate impressive performance benchmarks.  However, I found that the publication was a bit opaque with respect to how one might go about implementing it, and I found that most of the references cited were behind paywalls.  I put a good amount of my personal time into doing the legwork to implement their algorithm for uniform cubic B-splines in python (3.7 with numpy) and I would like to add my implementation to the body of literature available to speed things up for any curious readers like myself.
