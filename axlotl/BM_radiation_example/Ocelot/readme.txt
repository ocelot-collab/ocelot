The calculations were prepared for area x = [-10 mm, 10 mm] and y = [-10 mm, 10 mm].

bmrad.py - this module was changed by me for two-dimensional field calculation.
Ocelot_field.py - directly the code for the field calculation using formulas determined in bmrad.py

Ocelot_Field_Y_Npoints.txt or  -  horizontal and vertical polarization of calculated fields.
Ocelot_Field_X_Npoints.txt        Fields are presented in complex form.
                                  N in the names of the fields is a number of electron positions and horizontal angular deflections for which the fields were calculated.
								                  It means that in fact the fields were calculated for N^3 points.
								                  Vertical angular deflection was not taken into account because the calculations were done for the PETRA.
								                  The case when 1 point is used means calculation for 1 electron.

IntX_Npoints.png and IntY_Npoints.png - Saved pictures for the case with N points.

IntX_One_Electron.png and IntY_One_Electron.png - Saved pictures for one electron.
