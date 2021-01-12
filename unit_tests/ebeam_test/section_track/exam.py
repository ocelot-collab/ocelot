import sys
sys.path.append("/Users/tomins/Nextcloud/DESY/repository/ocelot/")
from ocelot import *
from accelerator.s2e_sections.sections import *
from ocelot.utils.section_track import *

data_dir = "unit_tests/ebeam_test/section_track/data"
tws0 = Twiss()
tws0.E = 0.005
tws0.beta_x = 0.286527307369
tws0.beta_y = 0.286527307369
tws0.alpha_x = -0.838833736086
tws0.alpha_y = -0.838833736086

all_sections = [A1, AH1, LH, DL, BC0]  # , L1, BC1]
parameter = 1
section_lat = SectionLattice(sequence=all_sections, tws0=tws0, data_dir=data_dir)
sec = all_sections[parameter]
s = section_lat.dict_sections[sec]
lat = s.lattice
if parameter == 0:
    energy = 0.005
elif parameter == 1:
    energy = 0.15
else:
    energy = 0.13

r_matrix = lattice_transfer_map(lat, energy=0.1)
print(r_matrix)

p_array = generate_parray()

section_lat.track_sections(sections=[A1], p_array=p_array)