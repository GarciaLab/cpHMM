Each file contains the source data for gap genes (Hb, Kr, Gt, Kni) processed dorsal gene expression levels measured using immunofluorescence techniques in 23 Drosophila embryos with FC depth comprised between 10 and 20 microns. We only picked embryos imaged close to their midsagittal plane (Orientation '1').

Column 1 contains the embryo number referring to its position on the slide (see Images.zip)
Column 2 contains an information about the azimuthal orientation of the embryo on the slide ('1' if the confocal plane is closer to midsagittal plane, '2' if it is closer to the coronal plane)
Column 3 contains the furrow canal depth (delta_FC) measured in micrometers.
Column 4 contains the corresponding estimated age in minutes.
Columns 5 to 1004 contain a 1x1000 vector representing the dorsal time-corrected gene expression level of the gap gene (Hb, Kr, Gt, Kni) in the embryo. The 1000 points are equally spaced along the AP axis. Thus, g345 represents the gene expression level at 34.5%EL. NaN means that we coudn't reliably detect the profile intensity at that position (usually near the edges).

Columns 1-4 are identical in the four files.