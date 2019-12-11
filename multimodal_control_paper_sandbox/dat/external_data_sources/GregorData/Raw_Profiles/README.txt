Each file contains the source data for raw gap genes (Hb, Kr, Gt, Kni) dorsal intensity profiles measured using immunofluorescence techniques in 163 Drosophila embryos during n.c. 14.

Column 1 contains the embryo number referring to its position on the slide (see Images.zip)
Column 2 contains an information about the azimuthal orientation of the embryo on the slide ('1' if the confocal plane is closer to midsagittal plane, '2' if it is closer to the coronal plane)
Column 3 contains the furrow canal depth (delta_FC) measured in micrometers.
Column 4 contains the corresponding estimated age in minutes.
Columns 5 to 1004 contain a 1x1000 vector representing the dorsal intensity profile of the gap gene (Hb, Kr, Gt, Kni) in the embryo. The 1000 points are equally spaced along the AP axis. Thus, Intensity345 represents the intensity at 34.5%EL. NaN means that we coudn't reliably detect the profile intensity at that position (usually near the edges).x

Columns 1-4 are identical in the four files.