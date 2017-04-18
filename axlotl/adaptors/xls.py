from openpyxl import Workbook
from ocelot import *
import ocelot

ignore_classes = (Edge,)

def create_xls(lat, file_name):
	
	wb = Workbook()

	ws = wb.get_active_sheet()


	for e in lat.sequence:
		
		if e.__class__ in ignore_classes:
			continue
		
		try:
			etype = e.type
		except:
			etype = None
			
		eclass = str(e.__class__).replace(str(e.__module__)+'.','').upper()
		
		ws.append([e.id, eclass, etype, e.l] )


	wb.save(file_name)


q1 = Quadrupole(l=0.1, k1=0.1, el_id = "Q1")
q1.type = "TYPEXXX"
q2 = Quadrupole(l=0.1, k1=-0.2, el_id = "Q2")
d = Drift(l=0.1)
b1 = SBend(l=0.1, angle = 0.001, el_id = "B1")

lat = MagneticLattice( (q1,d,q2,d,b1,d,b1))


create_xls(lat, file_name="test.xlsx")


