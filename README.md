# Multilateration
Draw loci corresponding to radio transmission multilateration.

This program plots loci for the scenario where there are radio towers
at known locations and a transmitter at an unknown location. The radio
towers accurately timestamp when they receive the transmission, allowing
time difference of arrival (TDOA) to be determined. This forms a
multilateration problem, producing n-1 loci where n is the number
of towers.

Only the 2-dimensional case is considered. It is assumed that the effect
on TDOA fron the vertical component of the transmission path is negligible.
For example, a path that is 5km horizontally and 500m vertically is
in total 5.025km ((5\**2 + 0.5\**2)\**0.5). Depending on clock noise this could
be considered negligible.
