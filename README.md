## The Spectral Analysis Routine for Asteroids, SARA

_Original Author: Sean S. Lindsay_

_E-mail:_ <Sean.Lindsay@physics.ox.ac.uk>
___
###General Information
The SARA code is written for the Interactive Data Language (IDL) v. 8.0+.  The code has been version tested back to IDL v.6.0 and to date shows no anomalous behavior on previous IDL versions.  The SARA code has also been tested for use on Windows, Mac, and Linux operating systems.  The main algorithm, SARA.pro, will run with no additional IDL packages installed.  The post-SARA analysis programs require the user to install the [IDL Astronomy Library](http:idlastro.gsfc.nasa.gov/) and the [IDL Coyote Library](http://www.idlcoyote.com/).

The SARA code performs a band parameter analysis for near infrared (NIR) reflectance spectral data for asteroids.  The band parameters are measured from the 1 and 2 &mu;m absorption bands due to crystal field absorptions from electronic transitions in Fe<sup>2+</sup> located in the M1 and/or M2 crystalography sites in olivines and pyroxenes.  Due to the lack of these absorption features in that majority of asteroid types, SARA is currently only configured to work for S- and V-type asteroids.  This version of SARA also includes four post-SARA analysis routines that make use of the band parameter analysis.  These programs need be run in sequence.  See below for details.

###SARA Package Contents

####sara_v01_users_guide.pdf
	The User's Manual for SARA.  This document contains information for how to set SARA up for use on your machine; provides referencing information; describes the algorithm methologies; and describes the functions of the four post-SARA analysis routines.

####SARA.pro
	The main SARA algorithm.  This program performs the 1 &mu;m and 2 &mu;m band parameter analysis and writes the band parameters, band definitions, band images, and an IDL save file to automatically created directories.  For initial use and keywords, please see the User's Guide.

####SARA_post1_master.pro
This post-SARA analysis program uses the SARA.pro outputs to collate a sample's band parameter individual band parameters into a pair of master files: one for all the polynomial orders, and one for the average of the three polynomial orders. For details, please see the User's Guide.

####SARA_post2_mineralogy.pro
This post-SARA analysis program uses the band parameter averaged master output file and performs mineralogical calculations using calibration equations from the following sources: S-type asteroids:  Dunn et al. (2010c); V-type asteroids: Burbine et al. (2007); generic type (S- and V-type): Gaffey et al. (2002)
	
####SARA_post3_analogs.pro
This post-SARA analysis program uses calculated band parameters to identify potential meteorite analogs to the asteroids.  This program generates a plot of the sample's Band I Centers versus Band Area Ratios (BAR) with known meteorite zones highlighted.  For details, please see the User's Guide.

####SARA_post4_ssubtypes.pro
This post-SARA analysis program uses calculated band parameters to identify the Gaffey et al. (1993) S-subtypes.  This program generates a plot of the sample's Band I Centers versus Band Area Ratios (BAR) overplotted on the seven [S(I)-S(VII)] S-subtype zones as well as the Basaltic Achondrite (BA) zone that is associated with HED meteorites.  For details, please see the User's Guide.

