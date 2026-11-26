import sys
W, H = 900, 380
INK, MUTED, CARD, BOX, BOXS, PANEL = '#1f2328', '#59636e', '#ffffff', '#e6eaef', '#c3cad3', '#f6f8fa'
REF, A, B, HI, HIBG = '#8c959f', '#cf222e', '#0969da', '#9a6700', '#fff1d6'
out = []

def t(x, y, s, size=16, fill=INK, w='600', anchor='middle'):
    out.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" font-weight="{w}" text-anchor="{anchor}">{s}</text>')

def rect(x, y, w, h, fill, stroke='none', rx=4, sw=1):
    out.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')

def path(d, stroke, sw=4):
    out.append(f'<path d="{d}" fill="none" stroke="{stroke}" stroke-width="{sw}" stroke-linecap="round" stroke-linejoin="round"/>')

def arrow(x1, x2, y):
    path(f'M{x1},{y} L{x2 - 4},{y}', MUTED, 2.5)
    out.append(f'<path d="M{x2},{y} l-12,-7 l0,14 Z" fill="{MUTED}"/>')

out.append(f'<rect width="{W}" height="{H}" rx="12" fill="{CARD}"/>')
t(W / 2, 52, 'graph-walks reads one stretch of a pangenome from a gfa-to-tabix walk file', 22)

PY, PH = 88, 240
LY = PY + PH + 26

P1, P1W = 30, 300
rect(P1, PY, P1W, PH, PANEL, BOXS, 12)
cx0, cx1, cy = P1 + 28, P1 + P1W - 28, 170
rect(cx0, cy - 12, cx1 - cx0, 24, BOX, BOXS, 12)
wx0, wx1 = cx0 + .58 * (cx1 - cx0), cx0 + .66 * (cx1 - cx0)
zx0, zx1, zy = cx0 + 20, cx1 - 20, 235
out.append(f'<path d="M{wx0:.1f},{cy + 12} L{zx0},{zy} L{zx1},{zy} L{wx1:.1f},{cy + 12} Z" fill="{HIBG}"/>')
rect(wx0, cy - 16, wx1 - wx0, 32, HI, rx=4)
rect(zx0, zy, zx1 - zx0, 36, HIBG, HI, 6, 2)
t((P1 + P1W / 2), 132, 'chromosome 6', 15, MUTED, '400')

P3, P3W = 420, 450
rect(P3, PY, P3W, PH, PANEL, BOXS, 12)
X0, X1, Y = P3 + 30, P3 + P3W - 30, 200
nodes = [(0, .2), (.24, .44), (.56, .76), (.8, 1)]
def nx(f): return X0 + f * (X1 - X0)
for a, b in nodes:
    rect(nx(a), Y - 16, nx(b) - nx(a), 60, BOX, BOXS, 4)
ins0, ins1 = nx(.45), nx(.55)
rect(ins0 - 6, Y - 70, ins1 - ins0 + 12, 30, '#ddeaff', B, 4, 1.5)
yb, yr, ya = Y, Y + 14, Y + 28
path(f'M{X0 - 8},{yb} L{nx(.44) - 14:.1f},{yb} C{nx(.44):.1f},{yb} {ins0 - 10:.1f},{Y - 55} {ins0 + 4:.1f},{Y - 55} L{ins1 - 4:.1f},{Y - 55} C{ins1 + 10:.1f},{Y - 55} {nx(.56):.1f},{yb} {nx(.56) + 14:.1f},{yb} L{X1 + 8},{yb}', B)
path(f'M{X0 - 8},{yr} L{X1 + 8},{yr}', REF)
path(f'M{X0 - 8},{ya} L{nx(.2) - 4:.1f},{ya} C{nx(.28):.1f},{ya + 62} {nx(.72):.1f},{ya + 62} {nx(.8) + 4:.1f},{ya} L{X1 + 8},{ya}', A)

arrow(P1 + P1W + 22, P3 - 22, PY + PH / 2)

t(P1 + P1W / 2, LY + 8, 'you pick a stretch', 17)
t(P3 + P3W / 2, LY + 8, 'you get each genome’s route', 17)

svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif">
<title>graph-walks: one stretch of a pangenome, read from a gfa-to-tabix walk file</title>
{chr(10).join(out)}
</svg>
'''
open(sys.argv[1], 'w').write(svg)
