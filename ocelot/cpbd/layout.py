import numpy as np
from ocelot.cpbd.magnetic_lattice import MagneticLattice
import plotly.graph_objects as go


class MachineLayout:
    def __init__(self):
        # List of beamlines with metadata
        # Structure: {'name': str, 'lattice': MagneticLattice, 'parent_name': str, 'anchor_element_id': str}
        self.lines = []
        self._surveys = {}  # Cache for calculated surveys {name: survey_list}

    def add_line(self, name, lattice, parent_name=None, anchor_element_id=None):
        """
        Adds a beamline to the facility.
        :param name: Unique name for this line (e.g. 'injector', 'dump_line')
        :param lattice: The MagneticLattice object
        :param parent_name: Name of the line this one branches FROM. None if it's the root.
        :param anchor_element_id: The ID of the element in the parent line AFTER which this line starts.
        """
        self.lines.append({
            'name': name,
            'lattice': lattice,
            'parent': parent_name,
            'anchor': anchor_element_id
        })

    def survey(self):
        """
        Computes surveys for all lines in parent-first order (independent of add_line order).
        """
        self._surveys = {}

        def get_angles_from_W(W):
            xpd, ypd, zpd = W[0, 2], W[1, 2], W[2, 2]
            phi = np.arcsin(np.clip(ypd, -1.0, 1.0))
            theta = np.arctan2(xpd, zpd)
            psi = np.arctan2(W[1, 0], W[1, 1])
            return theta, phi, psi

        # Index by name and validate uniqueness
        by_name = {}
        for line in self.lines:
            name = line["name"]
            if name in by_name:
                raise ValueError(f"Duplicate line name '{name}'")
            by_name[name] = line

        remaining = set(by_name.keys())

        # Resolve until all computed, or no progress
        while remaining:
            progressed = False

            for name in list(remaining):
                line = by_name[name]
                lat = line["lattice"]
                parent = line["parent"]
                anchor = line["anchor"]

                # Root line: compute immediately
                if parent is None:
                    mid, end = lat.survey(X0=0, Y0=0, Z0=0, theta0=0, phi0=0, psi0=0)
                    self._surveys[name] = end
                    remaining.remove(name)
                    progressed = True
                    continue

                # Non-root: parent must already be computed
                if parent not in self._surveys:
                    continue  # can't do this one yet

                parent_survey = self._surveys[parent]

                # Find anchor node in parent
                anchor_node = next(
                    (p for p in parent_survey if p["element"] and p["element"].id == anchor),
                    None
                )
                if anchor_node is None:
                    raise ValueError(f"Anchor '{anchor}' not found in parent line '{parent}' for child '{name}'")

                # Start from parent's anchor exit
                start_V = anchor_node["r_end"]
                start_W = anchor_node["W"]

                th0, ph0, ps0 = get_angles_from_W(start_W)

                mid, end = lat.survey(
                    X0=float(start_V[0]), Y0=float(start_V[1]), Z0=float(start_V[2]),
                    theta0=float(th0), phi0=float(ph0), psi0=float(ps0)
                )
                self._surveys[name] = end
                remaining.remove(name)
                progressed = True

            if not progressed:
                # We are stuck: cycle or missing parent
                missing = []
                for name in remaining:
                    parent = by_name[name]["parent"]
                    if parent not in by_name:
                        missing.append((name, parent))

                if missing:
                    msg = ", ".join([f"'{n}' -> missing parent '{p}'" for n, p in missing])
                    raise ValueError(f"Cannot compute survey: {msg}")

                # Otherwise it's a cycle (A->B->A) or chain with unresolved anchors
                deps = ", ".join([f"'{n}' depends on '{by_name[n]['parent']}'" for n in remaining])
                raise ValueError(f"Cannot compute survey due to cyclic dependencies or unresolved parents: {deps}")

        return self._surveys

    def check_interferences(self, min_distance=0.1):
        """
        Checks for collisions between different beamlines.
        :param min_distance: Minimum allowed distance [m] used if element width is 0.
        :return: List of collisions
        """
        collisions = []

        all_elements = []
        for line_name, survey_data in self._surveys.items():
            for item in survey_data:
                # Skip zero-length markers
                if item['element'] and item['element'].l > 0:
                    all_elements.append({
                        'line': line_name,
                        'el': item['element'],
                        'p0': item['r_start'],
                        'p1': item['r_end']
                    })

        import itertools
        for a, b in itertools.combinations(all_elements, 2):
            if a['line'] == b['line']:
                continue

            dist = self._segment_distance(a['p0'], a['p1'], b['p0'], b['p1'])

            # --- FIX STARTS HERE ---
            # 1. Access 'w', not 'width' (Your Element class stores 'w')
            # 2. Explicitly handle 0.0 width

            wa = getattr(a['el'], 'width', 0.0)

            if wa <= 0: wa = min_distance

            wb = getattr(b['el'], 'width', 0.0)
            if wb <= 0: wb = min_distance
            radius_sum = (wa + wb) / 2.0

            # Now: if dist is 0 and radius_sum is 0.1, 0 < 0.1 is True.
            if dist < radius_sum:
                collisions.append((a['line'], a['el'].id, b['line'], b['el'].id, dist))

        return collisions

    @staticmethod
    def _segment_distance(p1, p2, p3, p4):
        """
        Calculates minimum distance between two 3D segments (p1-p2) and (p3-p4).
        Standard algorithm.
        """
        u = p2 - p1
        v = p4 - p3
        w = p1 - p3

        a = np.dot(u, u)
        b = np.dot(u, v)
        c = np.dot(v, v)
        d = np.dot(u, w)
        e = np.dot(v, w)
        D = a * c - b * b

        sD = D
        tD = D

        SMALL_NUM = 1e-8

        if D < SMALL_NUM:  # Lines are parallel
            sN = 0.0
            sD = 1.0
            tN = e
            tD = c
        else:
            sN = (b * e - c * d)
            tN = (a * e - b * d)
            if sN < 0.0:
                sN = 0.0
                tN = e
                tD = c
            elif sN > sD:
                sN = sD
                tN = e + b
                tD = c

        if tN < 0.0:
            tN = 0.0
            if -d < 0.0:
                sN = 0.0
            elif -d > a:
                sN = sD
            else:
                sN = -d
            sD = a
        elif tN > tD:
            tN = tD
            if (-d + b) < 0.0:
                sN = 0
            elif (-d + b) > a:
                sN = sD
            else:
                sN = (-d + b)
            sD = a

        sc = 0.0 if abs(sN) < SMALL_NUM else sN / sD
        tc = 0.0 if abs(tN) < SMALL_NUM else tN / tD

        dP = w + (sc * u) - (tc * v)
        return np.linalg.norm(dP)




