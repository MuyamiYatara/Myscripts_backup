import numpy as np
import re
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType, Orbital
from scipy.interpolate import UnivariateSpline
from pymatgen.io.vasp.outputs import Outcar


class winzz(object):
    def __init__(self, file='vasprun.xml'):
        try:
            run = Vasprun(file)
            self.dos = run.tdos
            self.idos = run.idos
            self.energies = self.dos.energies
            self.pdos = run.pdos
            self.bandgap = run.eigenvalue_band_properties[0]
            self.atomic = (run.atomic_symbols)
            self.atomnum = len(self.atomic)
            self.atomtyp = len(set(self.atomic))
            self.atomcoord = run.structures[0].cart_coords
            self.projstr = ''
            self.winstr = ''
            self.nband = run.parameters['NBANDS']
            self.soc = run.parameters['LSORBIT']
            self.efermi = run.efermi
            self.atommap = list(zip(*[self.atomic, range(self.atomnum)]))

        except:
            print('\x1b[0;31;40m' +
                  'no VASPRUN file or VASPRUN is not correct!' + '\x1b[0m')

    def integrate(self, orb, start, end):
        x = self.energies
        y = orb
        g = UnivariateSpline(x, y, s=0)
        res = g.integral(start, end)
        return res

    def nume(self, enstart, enend):
        x = self.energies
        y = self.idos.densities[Spin.up]
        g = UnivariateSpline(x, y, s=0)
        if self.soc:
            return g(enend) - g(enstart)
        else:
            return (g(enend) - g(enstart)) / 2

    def output(self, tol=1, wtol=[0.4, 1.2, 2.0], froz=[None, None],  disen=[True, True], ex=[None, None]):
        def partition2(x):
            part = []
            leng = len(x) // 2
            for i in range(leng):
                part.append([x[2 * i], x[2 * i + 1]])
            return part

        extt = [(self.bandgap+.2)/2, (self.bandgap+.2)/2]

        if ex == [None, None]:
            ex = [0, 0]

        def ext(x): return [x[0] - extt[0] + ex[0], x[1] + extt[0] + ex[1]]

        def sub(x): return x[1] - x[0]

        dos = self.dos.densities[Spin.up]
        idos = self.idos.densities[Spin.up]
        engbut = min(self.energies)  # Minimum energy
        engtop = max(self.energies)  # Maximum energy
        idostop = max(idos)  # Maximum idos
        idosbut = min(idos)  # Minimum idos
        boundary = []
        iboundary = []

        # Energies for zero DOS

        for i in range(1, len(self.energies) - 1):
            if dos[i] == 0 and dos[i - 1] != 0:
                boundary.append(self.energies[i])
                iboundary.append(idos[i])
            if dos[i] == 0 and dos[i + 1] != 0:
                boundary.append(self.energies[i])
                iboundary.append(idos[i])

        if len(boundary) % 2 != 0:
            if dos[0] != 0:
                boundary = [engbut] + boundary
                iboundary = [idosbut] + iboundary
            else:
                boundary = boundary + [engtop]
                iboundary = iboundary + [idostop]

        boundary = list(map(ext, partition2(boundary)))
        iboundary = list(map(sub, partition2(iboundary)))
        bdict = dict(zip(iboundary, boundary))
        tmp = len(boundary)

        a = boundary
        ia = iboundary
        b = []
        ib = []
        ii = 0
        for begin, end in sorted(a):
            if b and b[-1][1] >= begin:
                b[-1] = (b[-1][0], end)
                ib[-1] = ib[-1] + ia[ii]
            else:
                b.append((begin, end))
                ib.append((ia[ii]))
            ii = ii + 1

        (boundary, iboundary) = (b, ib)

        bdict = dict(zip(iboundary, boundary))
        elenum = 0

        engwin = []
        for i in bdict:
            if bdict[i][0] <= self.efermi <= bdict[i][1]:
                engwin = bdict[i]
                elenum = i
                break

        if engwin == []:
            engwin = [0, 0]

            rebdict = dict(zip(list(reversed(iboundary)),
                               list(reversed(boundary))))

            for i in bdict:
                if bdict[i][0] > self.efermi:
                    engwin[1] = bdict[i][1]
                    elenum = i
                    break
            for i in rebdict:
                if rebdict[i][1] < self.efermi:
                    engwin[0] = bdict[i][0]
                    elenum = elenum + i
                    break

        engmin = engwin[0]
        engmax = engwin[1]

        def iondos(ion, orbital):
            if orbital == 's' or orbital == 0:
                orb = self.pdos[ion][Orbital.s][Spin.up]
            elif orbital == 'p' or orbital == 1:
                orb = self.pdos[ion][Orbital.px][Spin.up] + \
                    self.pdos[ion][Orbital.py][Spin.up] + \
                    self.pdos[ion][Orbital.pz][Spin.up]
            elif orbital == 'd' or orbital == 2:
                orb = self.pdos[ion][Orbital.dxy][Spin.up]+self.pdos[ion][Orbital.dyz][Spin.up] + \
                    self.pdos[ion][Orbital.dz2][Spin.up] + \
                    self.pdos[ion][Orbital.dxz][Spin.up] + \
                    self.pdos[ion][Orbital.dx2][Spin.up]
            return orb

        def wight(ion, orbital):
            orb = iondos(ion, orbital)
            return self.integrate(orb, engmin, engmax)

        def cstr(x): return str(x[0]) + ',' + str(x[1]) + ',' + str(x[2])

        nwann = 0
        proj = {}
        for i in range(self.atomnum):
            proj[i] = ''

        for i in range(self.atomnum):
            for j in range(3):
                print(i, j, wight(i, j))
                if wight(i, j) >= wtol[j]:
                    proj[i] = proj[i] + 'l=' + str(j) + ';'
                    if self.soc == True:
                        nwann += 2 * (2 * j + 1)
                    else:
                        nwann += (2 * j + 1)

        newtol = np.array(wtol)
        if not self.soc:
            elenum = elenum / 2

        while nwann > elenum:
            nwann = 0
            proj = {}
            newtol = newtol + tol * np.array([.1, .3, .5])
            print('\x1b[6;31;47m' + 'WARNING! INITITAL WTOL IS TOO SMALL SO THE BAND IS NOT ENOUGH TO CONSTRACT THE WANNIER FUCTIONS I WILL CHOSE LARGER WTOL, THE NEW TOL IS ' + str(newtol) + '\x1b[0m')
            for i in range(self.atomnum):
                proj[i] = ''
            for i in range(self.atomnum):
                for j in range(3):
                    if wight(i, j) > newtol[j]:
                        proj[i] = proj[i] + 'l=' + str(j) + ';'
                        if self.soc == True:
                            nwann += 2 * (2 * j + 1)
                        else:
                            nwann += (2 * j + 1)
        mengmax = engmax
        melenum = elenum

        while 1.2 * nwann < melenum:
            mengmax = mengmax - 0.1
            melenum = self.nume(engmin, mengmax)
            if mengmax < self.efermi:
                break

        engmax = mengmax

        for i in list(proj):
            if proj[i] != '':
                proj[i] = 'c=' + cstr(self.atomcoord[i]) + ':' + proj[i][:-1] + \
                    '   !  ' + str(i + 1) + '  ' + self.atomic[i] + '\n'
            else:
                proj.pop(i, None)
        projstr = 'Begin Projections\n'
        for i in proj:
            projstr += proj[i]
        projstr += 'End Projections\n'
        projstr += 'num_wann=' + str(nwann)
        projstr += '\nnum_bands=' + str(self.nband)
        dismin = disen[0]
        dismax = disen[1]
        if dismin:
            projstr += '\ndis_win_min=' + str(engmin)
        if dismax:
            projstr += '\ndis_win_max=' + str(engmax)
        projstr += '''
num_iter = 100
dis_num_iter = 500
search_shells = 25
iprint = 2
dis_conv_tol = 1.0E-10
num_cg_steps = 30
guiding_centres = .true.'''
        if self.soc:
            projstr += '\nspinors=.True.'

        frozmin = froz[0]
        frozmax = froz[1]
        projstr += '\nfermi_energy =' + str(self.efermi)

        if frozmin == True:
            projstr += '\ndis_froz_min=' + str(engmin)
        elif frozmin is not None:
            projstr += '\ndis_froz_min=' + str(self.efermi - frozmin)

        if frozmax is not None:
            projstr += '\ndis_froz_max=' + str(self.efermi + frozmax)

        with open('wannier90.win_auto', 'w') as f:
            f.write(projstr)

        if nwann == 0:
            print('Error')


def convert_kpoints_format(kpoints_file, win_file):
    pattern = re.compile(
        r'^\s*([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+!\s*(.+)')

    with open(kpoints_file, 'r') as f:
        lines = f.readlines()

    points = []
    for line in lines:
        m = pattern.match(line)
        if m:
            x, y, z = m.group(1), m.group(2), m.group(3)
            label = m.group(4).strip()
            label = label.lstrip('\\')
            points.append((label, x, y, z))

    if not points:
        print("No matching coordinate lines found. Check your KPOINTS format!")
        return

    segments = []
    for i in range(0, len(points) - 1, 2):
        segments.append((points[i], points[i+1]))

    output_lines = [
        "bands_plot=.true.",
        "bands_num_points=30",
        "begin kpoint_path"
    ]
    for (lbl1, x1, y1, z1), (lbl2, x2, y2, z2) in segments:
        output_lines.append(f"{lbl1} {x1} {y1} {z1} {lbl2} {x2} {y2} {z2}")
    output_lines.append("end kpoint_path")

    with open(win_file, 'a') as f:
        f.write("\n" + "\n".join(output_lines) + "\n")

    print(f"Converted {len(segments)} segments and appended to {win_file}")


if __name__ == '__main__':
    pass

a = winzz('./vasprun.xml')
a.output(froz=[True, 1])
# convert_kpoints_format("../band/KPOINTS", "wannier90.win")
