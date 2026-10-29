# Mn3Sn annotated notebook: digital transcription

## Status and use

This is a page-addressable digital transcription of [Mn3Sn.pdf](../Mn3Sn.pdf),
with the same local notebook format as the Alpha-Fe2O3 extraction. Open the
[interactive notebook](Mn3Sn-notebook.html) to compare the text with the original
color page images and record annotations.

The 23 retained notebook pages are represented in order. Equations, numerical
tables, readable notes, and diagram labels are transcribed below; drawings and
pasted figures are described and retained in full in the page images. This is
a semantic transcription, not a pixel-by-pixel redraw or independently verified
model specification. The source PDF remains authoritative for handwriting and
graphical detail. See [README.md](../README.md) and [plan.md](../plan.md) for the
separate description and audit of the MATLAB implementation.

Conventions:

- Page numbers are PDF page numbers, starting at 1. IDs such as
  `P07-DMI-FIELDS` provide stable discussion targets.
- Quoted wording and equations record the notebook, including unresolved
  questions, omitted factors, incomplete steps, and inconsistent conventions.
  No formula has been silently repaired or adopted as a simulation definition.
- `[unclear]` marks a reading that cannot be established confidently.
  `UNRESOLVED` marks an interpretation or internal discrepancy left open.
- Colors are local diagram notation. Their meaning can change between pages;
  there is no global mapping from ink color to chemical species or sublattice.
- Source labels `[1]`–`[5]` refer to the notebook's own abbreviated reference
  list on PDF page 23, not the numbered publications in the repository README.
- Axis directions, state angles, and labels are local to each page. Coordinate
  transformations are recorded as drawn, without merging different conventions.
- In the array transcriptions, rows run top to bottom and columns left to right.
  `X` records a drawn cross; `0`, `1`, `3`, and `A` retain the notebook's symbols.
  A cross is not silently replaced by a numerical zero.
- Author explanations and code correspondence are initially unfilled. The HTML
  notebook provides separate annotation fields, saved in browser storage when
  available, with a JSON export. Annotations do not edit this Markdown file.

Extraction record: 2026-10-21; Poppler `pdftoppm`, 160 dpi, JPEG quality 90;
24 full-color page images. PDF pages 19 and 22 are landscape and retain their
original orientation. Page 4's image and transcription were refreshed on
2026-10-21 for the added spherical-coordinate formula and angle annotations.
Re-rendering all 24 pages with the command below confirmed that the other
23 JPEGs are byte-identical to the original extraction.
Original extraction source size: 25,963,296 bytes. SHA-256:
`39bc60d99fc85b32f7b44870a20aa08c7121f5da962d0d802366952c0902a1d1`.

Re-rendering command, run from the repository root after creating the image
directory (overwrites the generated page JPEGs):

```bash
pdftoppm -r 160 -jpeg -jpegopt quality=90 Mn3Sn.pdf docs/Mn3Sn-pages/page
```

The PDF text layer was inspected but is unreliable for handwritten mathematics.
Every page was visually reviewed. This extraction does not verify the notebook's
physics, identify incomplete references, or compare its masks against MATLAB.

Page removal: 2026-10-27. Original page 23 was removed from the HTML notebook
and this transcript; References is now page 23. Its stable ID `P24-REFERENCES`,
image filename `page-24.jpg`, and browser annotation key `P24` retain their
original numbering to preserve existing references and saved notes. PDF links
use the revised page numbering.

## Page index

| PDF page | Stable ID | Notebook section |
|---:|---|---|
| 1 | `P01-CRYSTAL-STRUCTURE` | Structure, direction indices, easy-plane questions |
| 2 | `P02-STABLE-STATES` | Six configurations and source comparisons |
| 3 | `P03-LATTICE-SPINS` | Array region, spin pattern, local easy-axis sketch |
| 4 | `P04-ANISOTROPY` | Anisotropy derivation and Configuration II axes |
| 5 | `P05-ANISOTROPY-OLD` | Older rotated-axis derivation and field signs |
| 6 | `P06-DMI-ENERGY` | Cyclic three-sublattice energy |
| 7 | `P07-DMI-FIELDS` | DMI fields with the D vector along x |
| 8 | `P08-SOT-GEOMETRIES` | Configurations (a), (b), and sub-configuration question |
| 9 | `P09-SOT-GEOMETRIES` | Configurations (b), (c), with and without field |
| 10 | `P10-SWITCHING-POLARITY` | Current/field signs and six state angles |
| 11 | `P11-AXIS-TRANSFORM` | Axis transformation and transformed switching table |
| 12 | `P12-SWITCHING-SUMMARY` | Summary of the transformed sign table |
| 13 | `P13-CONFIGURATION-I` | Zhengde Configuration I state selection |
| 14 | `P14-AXIS-RELABELING` | Further axis and angle relabeling |
| 15 | `P15-CONFIGURATION-II` | Six-by-four switching table |
| 16 | `P16-SIX-STATE-SKETCH` | Six configurations and drive directions |
| 17 | `P17-PARAMETERS-UNITS` | Field-response sketches, parameters, unit checks |
| 18 | `P18-IN-PLANE-FIELD` | Configuration I sketches under opposite fields |
| 19 | `P19-ATOMISTIC-5X6` | Version 3 nearest-neighbor sketch, 5 × 6 example |
| 20 | `P20-OCCUPATION` | Lattice drawing and atom-type array |
| 21 | `P21-ANISOTROPY-EXCHANGE-MASKS` | Anisotropy and eight exchange masks |
| 22 | `P22-ATOMISTIC-7X6` | Version 3 nearest-neighbor sketch, 7 × 6 example |
| 23 | `P24-REFERENCES` | Five abbreviated references |

## Page 1 — Crystal structure and direction indices

Source: [PDF page 1](../Mn3Sn.pdf#page=1) · [page image](Mn3Sn-pages/page-01.jpg)

### `P01-CRYSTAL-STRUCTURE`

The heading is “crystal structure.” The upper figure is labeled “without HM”
and shows a kagome-like arrangement with red and gray sites, blue spin arrows,
bond labels 2.83 and 2.84, and printed crystallographic directions. Green
overlays introduce axes $a$, $b$, $a_2$, and $a_3$. No unit is visible beside
the two bond-length labels in this clipped figure.

Red questions:

> 没有考虑 easy plane anisotropy？导致 simulation 电流过大？

Translation: “Was easy-plane anisotropy not considered? Does this cause the
simulation current to be too large?”

> does the theory of reported anisotropy wrong? How to introduce easy-plane anisotropy?

Status: `UNRESOLVED`; these are questions in the source, not conclusions.

### `P01-DIRECTION-CONVERSION`

The notebook writes the conversion $[uvw]\Rightarrow[UVTW]$ as

$$
U=\frac{2u-v}{3},\qquad V=\frac{2v-u}{3},\qquad
T=-(U+V),\qquad W=w.
$$

The first worked tuple is

$$
(1/2,-1,0)\Rightarrow[2/3,-5/6,1/6,0].
$$

Below, axis sketches label $a$, $b$, $c$, and $d$. An orange note reads
“全部用 $\vec a,\vec b$ 来表示” (“express everything using $\vec a,\vec b$”).
The subsequent written chains are:

$$
c:\quad[-1/2,-1/2,0]
\Rightarrow[-1/3,-1/3,2/3,0]
\Rightarrow[2,-1,-1,0],
$$

$$
a:\quad[1,0,0]\Rightarrow[2,-1,-1,0]
\Rightarrow[-1,2,-1,0],
$$

$$
d:\quad[-1,2,0]\Rightarrow[0,1,-1,0]
\Rightarrow[-1,0,1,0].
$$

The final arrows are retained literally; their relation to axis relabeling or
symmetry is `UNRESOLVED`. They should not be read as verified numerical
equalities under a single unchanged basis.

### `P01-WITH-HM`

The lower section is labeled “with HM [4].” A pasted Mn3Sn/W/MgO device figure
contains current, field, strain, spin, and crystallographic arrows. The clipped
paragraph describes a hexagonal $D0_{19}$ structure and gives the orientation
relationship

$$
\mathrm{MgO}(110)[001]\parallel
\mathrm{Mn}_3\mathrm{Sn}(01\bar10)[0001].
$$

The paragraph also mentions cross-sectional TEM. It is a source excerpt in
the notebook; the referenced article was not independently checked here.

## Page 2 — Stable states

Source: [PDF page 2](../Mn3Sn.pdf#page=2) · [page image](Mn3Sn-pages/page-02.jpg)

### `P02-STABLE-STATES`

Six hand-drawn hexagonal motifs compare black and orange triangles with green
spin arrows. The source-[4] circled labels appear in this order:

| Position | Left | Middle | Right |
|---|---|---|---|
| Top row | 1 | 2 | 4 |
| Bottom row | 6 | 5 | 3 |

Below, a source-[1] figure displays the pair $+T_x^\gamma$ and $-T_x^\gamma$.
A source-[4] panel shows six states labeled 1–6 with $\varepsilon=0$ and
crystallographic direction arrows. Hand sketches are labeled “[2] Fig 1(d)”
and “[3].” The lower part begins a lattice drawing with alternating red and
black outlined circles on slanting lines.

The drawings are preserved in the image. No common angular origin or
one-to-one correspondence between all four source conventions is specified
on this page.

## Page 3 — Lattice region and spin pattern

Source: [PDF page 3](../Mn3Sn.pdf#page=3) · [page image](Mn3Sn-pages/page-03.jpg)

### `P03-LATTICE-SPINS`

The upper drawing extends the red/black circle lattice and outlines a
rectangular region. The middle drawing overlays horizontal and diagonal
lines and a repeating blue three-direction spin pattern. The three arrow
directions are also drawn separately at the lower right of this pattern.

At the bottom, a source-[3] figure labeled “Real space” includes directions
$[01\bar10]$, $[2\bar1\bar10]$, and $[0001]$. Beside it, a triangular sketch
labels Mn, Sn, and orange $H_k$ arrows. Blue/green arrows indicate local spins.
The source image retains their exact positions.

Status: the rectangular selection and color classes are graphical notebook
information; exact metric coordinates and correspondence to code indices
are not supplied here.

## Page 4 — Anisotropy

Source: [PDF page 4](../Mn3Sn.pdf#page=4) · [page image](Mn3Sn-pages/page-04.jpg)

### `P04-ANISOTROPY`

Three sketches label spins $\vec m_a,\vec m_b,\vec m_c$ and local axes
$\vec k_a,\vec k_b,\vec k_c$. The written energy and final derivative are

$$
E=k_u\sin^2\theta=k_u(1-\cos^2\theta)
=k_u-k_u(\vec k_a\cdot\vec m)^2,
$$

$$
\vec H_a=-\frac{\partial E}{\partial\vec m}
=2k_u(\vec k_a\cdot\vec m_a)\vec k_a.
$$

The intermediate handwritten fraction has $k_u(\vec k_a\cdot\vec m)^2$
over $\partial\vec m$; its numerator derivative symbol is omitted. “Similarly”:

$$
\vec H_b=2k_u(\vec k_b\cdot\vec m_b)\vec k_b,\qquad
\vec H_c=2k_u(\vec k_c\cdot\vec m_c)\vec k_c.
$$

No moment or magnetization denominator appears in these lines. This is the
notebook's notation, not a declaration that these expressions have tesla units.

### `P04-CONFIGURATION-II-AXES`

The lower heading is “zhengde paper / Configuration II.” The axis diagram has
$z$ upward, $y$ rightward, and $x$ out of the page; $P$ is drawn out of the
page and $I\parallel H$ rightward. The state is labeled $f$, $330^\circ$.
Above the easy-axis list, the revised page writes the Cartesian components

$$
[\sin\theta\cos\varphi,\sin\theta\sin\varphi,\cos\theta].
$$

The listed easy axes are

$$
\vec k_1=\vec k_a=[0,\sin(\pi/6),\cos(\pi/6)],
$$

$$
\vec k_2=\vec k_b=[0,\sin(3\pi/2),\cos(3\pi/2)],
$$

$$
\vec k_3=\vec k_c=[0,\sin(5\pi/6),\cos(5\pi/6)].
$$

The added notes beside the three vectors specify:

| Easy axis | $\theta$ | $\varphi$ |
|---|---|---|
| $\vec k_1=\vec k_a$ | $\pi/6$ | $\pi/2$ |
| $\vec k_2=\vec k_b$ | $3\pi/2$ | $\pi/2$ |
| $\vec k_3=\vec k_c$ | $5\pi/6$ | $\pi/2$ |

Reading these components with $\varphi=\pi/2$ gives
$[0,\sin\theta,\cos\theta]$ in the drawn $yz$ plane. The notebook retains
$\theta=3\pi/2$ for $\vec k_b$, giving $[0,-1,0]$; this is its full-turn
in-plane angle convention. The $\theta$ in the upper energy sketch instead
labels the angle between the spin and its easy axis. The revision makes the
coordinate convention explicit without changing the three listed vectors.

## Page 5 — Anisotropy: old discussion

Source: [PDF page 5](../Mn3Sn.pdf#page=5) · [page image](Mn3Sn-pages/page-05.jpg)

### `P05-ANISOTROPY-OLD`

The heading explicitly says “anisotropy — old discussion.” The first sketches
show a spin at angle $\theta$ to $\hat z$, a rotated axis $\hat z'$, and fields
$2k_u\cos\theta\,\hat z$ and $2k_u\cos\theta'\,\hat z'$.
The rotated-axis expansion reads

$$
2k_u\cos(45^\circ-\theta)\frac{\sqrt2}{2}(\hat z+\hat x)
$$

$$
=\sqrt2k_u(\cos45^\circ\cos\theta+
\sin45^\circ\sin\theta)(\hat z+\hat x)
$$

$$
=(k_u\cos\theta+k_u\sin\theta)(\hat z+\hat x).
$$

The lower derivation writes a positive derivative of a negative anisotropy
energy, with the resulting negative field:

$$
\vec H_a=
\frac{\partial[-k(\vec k_a\cdot\vec m_a)^2]}
{M_s\,\partial\vec m_a}
=-\frac{2k}{M_s}(\vec k_a\cdot\vec m_a)\vec k_a.
$$

A check mark is drawn beside it. At the bottom the page also writes

$$
1-\cos^2\theta=1-(\hat z\cdot\vec m)^2,
\qquad 2k\cos\theta\,\hat z,
$$

and, on the last line,

$$
\vec H_a=\frac{2k_u}{M_s}\cos\theta\,\hat z.
$$

Status: `UNRESOLVED`. The negative field line and positive final expression
are both retained; the sign and normalization are not reconciled here.

## Page 6 — DMI and three-sublattice energy

Source: [PDF page 6](../Mn3Sn.pdf#page=6) · [page image](Mn3Sn-pages/page-06.jpg)

### `P06-DMI-ENERGY`

The heading is “DMI,” with a triangle labeled $A$ at the top, $B$ lower right,
and $C$ lower left. Source label `[5]` precedes

$$
U=D\hat e_z\cdot\sum_{\langle\mu\nu\rangle}
\vec m_\mu\times\vec m_\nu.
$$

Literal note:

> $\langle\mu\nu\rangle$ indicates sum over $(\mu\nu)=(A,B),(B,C)$ and $(C,A)$.

The explicit expansion is

$$
U=D\hat e_z\cdot(\vec m_A\times\vec m_B+
\vec m_B\times\vec m_C+\vec m_C\times\vec m_A).
$$

### `P06-COMBINED-ENERGY`

Further down, the notebook writes

$$
U=J_0\sum_{\langle\mu\nu\rangle}\vec m_\mu\cdot\vec m_\nu
+D\hat e_z\cdot\sum\vec m_\mu\times\vec m_\nu
-k\sum_{\mu=A,B,C}(\vec m_\mu\cdot\vec e_{k,\mu})^2.
$$

Its next displayed lines expand only the exchange and DMI terms:

$$
J_0(\vec m_A\cdot\vec m_B+\vec m_B\cdot\vec m_C+
\vec m_C\cdot\vec m_A)
$$

$$
+D\hat e_z\cdot(\vec m_A\times\vec m_B+
\vec m_B\times\vec m_C+\vec m_C\times\vec m_A).
$$

The anisotropy term is not carried into that visible expansion. The omitted
continuation is not reconstructed.

## Page 7 — DMI fields for a vector along x

Source: [PDF page 7](../Mn3Sn.pdf#page=7) · [page image](Mn3Sn-pages/page-07.jpg)

### `P07-DMI-FIELDS`

Literal assumption: “Assume D along $\hat x$ direction.” The page changes
the perpendicular direction relative to page 6 and writes

$$
E=D\hat x\cdot(\vec m_1\times\vec m_2+
\vec m_2\times\vec m_3+\vec m_3\times\vec m_1).
$$

The stated scalar triple-product identity is

$$
\vec a\cdot(\vec b\times\vec c)
=\vec b\cdot(\vec c\times\vec a)
=\vec c\cdot(\vec a\times\vec b).
$$

The three regroupings isolate each spin:

$$
E=D[\vec m_1\cdot(\vec m_2\times\hat x)
+\hat x\cdot(\vec m_2\times\vec m_3)
+\vec m_1\cdot(\hat x\times\vec m_3)],
$$

$$
E=D[\vec m_2\cdot(\hat x\times\vec m_1)
+\vec m_2\cdot(\vec m_3\times\hat x)
+\hat x\cdot(\vec m_3\times\vec m_1)],
$$

$$
E=D[\hat x\cdot(\vec m_1\times\vec m_2)
+\vec m_3\cdot(\hat x\times\vec m_2)
+\vec m_3\cdot(\vec m_1\times\hat x)].
$$

With the notebook's definition $H_{\mathrm{DMI},i}=-\partial E/\partial m_i$:

$$
\vec H_{\mathrm{DMI},1}
=-D(\vec m_2\times\hat x+\hat x\times\vec m_3)
=-D[0,m_{2z}-m_{3z},-m_{2y}+m_{3y}],
$$

$$
\vec H_{\mathrm{DMI},2}
=-D(\hat x\times\vec m_1+\vec m_3\times\hat x)
=-D[0,m_{3z}-m_{1z},m_{1y}-m_{3y}],
$$

$$
\vec H_{\mathrm{DMI},3}
=-D(\hat x\times\vec m_2+\vec m_1\times\hat x)
=-D[0,m_{1z}-m_{2z},m_{2y}-m_{1y}].
$$

Two determinant calculations show
$\vec m_2\times\hat x=\hat y m_{2z}-\hat z m_{2y}$ and
$\hat x\times\vec m_3=-\hat y m_{3z}+\hat z m_{3y}$.
One intermediate coordinate-expansion line under the first field omits the
$-D$ prefactor, which is present again in the following lines. No moment or
magnetization divisor is written on this page.

## Page 8 — SOT switching geometries (a) and (b)

Source: [PDF page 8](../Mn3Sn.pdf#page=8) · [page image](Mn3Sn-pages/page-08.jpg)

### `P08-SOT-GEOMETRIES`

The heading is “SOT switching.” A red question reads:

> 算黑的3个和 orange 的3个会有区别吗？（boyu project, rotate plane）

Translation: “Would calculating the three black sites and the three orange
sites give different results? (Boyu project, rotate plane).”

Each group contains six green spin configurations on black/orange triangles.

| Group on this page | Written geometry | Drive marks |
|---|---|---|
| Configuration (a) in [3], top | $I\perp$ kagome plane; $\sigma\parallel$ kagome plane | $I,H$ into the page; blue polarization rightward |
| Configuration (a) in [3], middle | Same written geometry | $I$ into the page; no $H$ label |
| Configuration (b) in [3], bottom | $I\parallel$ kagome plane; $\sigma\perp$ kagome plane | $I,H$ leftward; blue polarization into the page |

The upper axes are annotated “Following [4]”: $z$ upward along
$[01\bar10]$, $x$ rightward along $[2\bar1\bar10]$, and $y$ into the page
along $[0001]$. The bottom group draws a different axis triad.
The handwritten polarization symbol is rendered as $\sigma$ here; its
arrow/circle marks are retained in the source image.

Status: the black-versus-orange question is `UNRESOLVED`. These are geometry
sketches, not measured switching outcomes.

## Page 9 — SOT switching geometries (b) and (c)

Source: [PDF page 9](../Mn3Sn.pdf#page=9) · [page image](Mn3Sn-pages/page-09.jpg)

### `P09-SOT-GEOMETRIES`

Three further groups each contain six spin configurations:

| Group on this page | Written geometry | Drive marks |
|---|---|---|
| Configuration (b) in [3], top | $I\parallel$ kagome plane; $\sigma\perp$ kagome plane | $I$ leftward; no $H$ label |
| Configuration (c) in [3], middle | $I\parallel$ kagome plane; $\sigma\parallel$ kagome plane | $I,H$ upward; blue polarization rightward |
| Configuration (c) in [3], bottom | Same written geometry | $I$ upward; no $H$ label |

The middle-group axes have $x$ rightward, $y$ upward, and $z$ out of the page.
Together, pages 8–9 record field-labeled and current-only versions of the three
geometries; no time evolution or numerical switching probability is supplied.

## Page 10 — Switching polarity and state angles

Source: [PDF page 10](../Mn3Sn.pdf#page=10) · [page image](Mn3Sn-pages/page-10.jpg)

### `P10-SWITCHING-POLARITY`

At the top:

$$
J_s=\theta\,\sigma\times J_c.
$$

The direction labels underneath are $J_s:-z$, $\sigma:+y$, $J_c:+x$.
Source label `[3]` accompanies a spin triangle, $\theta_{\mathrm{SH}}>0$,
$I,H\parallel x$ into the page, and $\vec P$ leftward. The upper example is
marked $330^\circ$; another is marked $30^\circ$. Handwritten cross-references
say “Corresponding to ③ in page 1” and “Corresponding to ⑤ in page 1.” These
are retained as historical references; they are not resolved to PDF page 1.
A red “why?” appears below the polarization sketch.

The current/field table drawn at the upper right is:

| Current sign | Field sign | Written final angle |
|---|---|---|
| $+J$ | $H_x>0$ | $330^\circ$ |
| $+J$ | $H_x<0$ | $30^\circ$ |
| $-J$ | $+H_x$ | $30^\circ$ |
| $-J$ | $-H_x$ | $330^\circ$ |

### `P10-SIX-ANGLES`

The six spin-triangle sketches below have this spatial order:

| Position | Left | Middle | Right |
|---|---|---|---|
| Top row | $330^\circ$ | $150^\circ$ | $210^\circ$ |
| Bottom row | $30^\circ$ | $90^\circ$ | $270^\circ$ |

Two pasted Mn3Sn/W device illustrations are captioned “Switching the same as
STT” and “type x switching.” These are notebook comparisons; their physical
equivalence is not established by this transcription.

## Page 11 — Transforming the axes

Source: [PDF page 11](../Mn3Sn.pdf#page=11) · [page image](Mn3Sn-pages/page-11.jpg)

### `P11-AXIS-TRANSFORM`

Literal note:

> To compare with zhengde results of Conf. I, we need to transform the axis

The three stages retain $J_s=\theta\,\sigma\times J_c$ and write these
direction triples:

| Stage | $J_s$ | $\sigma$ | $J_c$ |
|---|---|---|---|
| Left | $-z$ | $+y$ | $+x$ |
| Middle | $y$ | $-z$ | $-x$ |
| Right | $+y$ | $-z$ | $+x$ |

Orange note: “坐标系改变，$\Rightarrow x$ 反号” (“the coordinate system changes,
so the sign of x reverses”). Blue note: “$-y$ 变 $+y$，可以看成 $J_c$ 反号”
(“changing −y to +y can be regarded as reversing $J_c$”). In the rightmost
sketch, $x$ points out of the page, $y$ down, and $z$ right. The spin labels
are $m_1$ upper left, $m_2$ right, $m_3$ lower left. It is marked $210^\circ$,
with $H\parallel-x$, $I$ along $+x$, and $\vec P$ leftward.

The middle six-state drawings give:

| Position | Left | Middle | Right |
|---|---|---|---|
| Top row | $210^\circ$, d | $30^\circ$, a | $330^\circ$, f |
| Bottom row | $150^\circ$, c | $90^\circ$, b | $270^\circ$, e |

### `P11-TRANSFORMED-SWITCHING`

At the bottom, the page-10 table is redrawn at left. An arrow is annotated
“I 不变 / H 反号” (“I unchanged / H changes sign”). The right-hand table reads:

| Current | Field | Written final angle |
|---|---|---|
| $J_x>0$ | $H_x<0$ | $210^\circ$ |
| $J_x>0$ | $H_x>0$ | $150^\circ$ |
| $J_x<0$ | $H_x<0$ | $150^\circ$ |
| $J_x<0$ | $H_x>0$ | $210^\circ$ |

The page combines axis relabeling and current/field sign reasoning. No unique
rotation matrix is supplied; one is not inferred here.

## Page 12 — Switching summary

Source: [PDF page 12](../Mn3Sn.pdf#page=12) · [page image](Mn3Sn-pages/page-12.jpg)

### `P12-SWITCHING-SUMMARY`

The heading is “Summary.” The left sketch repeats the final convention from
page 11: $J_s=\theta\,\sigma\times J_c$ with directions $(+y,-z,+x)$,
$\theta_{\mathrm{SH}}>0$, $H\parallel-x$, $I$ along $+x$ out of the page,
and $\vec P$ leftward. The starting sketch is marked $210^\circ$.

| Current | Field | Written final angle |
|---|---|---|
| $J_x>0$ | $H_x<0$ | $210^\circ$ |
| $J_x>0$ | $H_x>0$ | $150^\circ$ |
| $J_x<0$ | $H_x<0$ | $150^\circ$ |
| $J_x<0$ | $H_x>0$ | $210^\circ$ |

This page repeats the notebook's sign summary. It supplies no numerical trace
or simulation settings to validate those outcomes.

## Page 13 — Zhengde paper: Configuration I

Source: [PDF page 13](../Mn3Sn.pdf#page=13) · [page image](Mn3Sn-pages/page-13.jpg)

### `P13-CONFIGURATION-I`

The axes have $z$ up, $x$ right, and $y$ into the page. The relation
$J_s=\theta\,\sigma\times J_c$ has direction labels $(z,x,y)$ underneath.
$P$ points right and $I\parallel H$ is drawn into the page. The six states
are ordered:

| Position | Left | Middle | Right |
|---|---|---|---|
| Top row | f, $330^\circ$ | a, $30^\circ$ | b, $90^\circ$ |
| Bottom row | c, $150^\circ$ | d, $210^\circ$ | e, $270^\circ$ |

Two initial states, a ($30^\circ$) and f ($330^\circ$), are followed by the
same two outcome sketches:

| Initial state | Written current/field combinations | Written outcome |
|---|---|---|
| a or f | $-J_c,+H_{\mathrm{ext}}$; $+J_c,-H_{\mathrm{ext}}$ | a, $30^\circ$ |
| a or f | $-J_c,-H_{\mathrm{ext}}$; $+J_c,+H_{\mathrm{ext}}$ | f, $330^\circ$ |

No outcome rows for initial states b–e are drawn on this page.

## Page 14 — Further axis and angle relabeling

Source: [PDF page 14](../Mn3Sn.pdf#page=14) · [page image](Mn3Sn-pages/page-14.jpg)

### `P14-AXIS-RELABELING`

Literal note:

> To unify the axis, we need to transform the axis

The relation $J_s=\theta\,\sigma\times J_c$ is annotated with $(+z,-y,+x)$.
The triad shows $z$ up, $-y$ left, and $x$ into the page. The top spin sketch
labels $m_2$ at the top, $m_3$ lower left, and $m_1$ lower right.

This sparse page records angle labels at positions corresponding to a
six-state layout:

| Position | Left | Middle | Right |
|---|---|---|---|
| Top row | $30^\circ$ | $330^\circ$ | $270^\circ$ |
| Bottom row | $210^\circ$ | $150^\circ$ | $90^\circ$ |

Two additional lines read, left to right,
$330^\circ,330^\circ,30^\circ$ and $30^\circ,330^\circ,30^\circ$.
Their inputs and transition arrows are not explicitly labeled on this page.
Status: `UNRESOLVED`; no missing switching table is reconstructed.

## Page 15 — Configuration II switching table

Source: [PDF page 15](../Mn3Sn.pdf#page=15) · [page image](Mn3Sn-pages/page-15.jpg)

### `P15-CONFIGURATION-II`

The triad has $z$ upward, $y$ rightward, and $x$ out of the page.
$P$ points out of the page and $I\parallel H$ rightward. The six states
use a = $30^\circ$, b = $90^\circ$, c = $150^\circ$, d = $210^\circ$,
e = $270^\circ$, and f = $330^\circ$. On the f sketch the spin labels are
$m_1$ lower right, $m_2$ top, and $m_3$ lower left.

The complete table is transcribed below. The four columns preserve the
notebook's order, including the same outcomes in all six rows:

| Initial state | $J_{+y},H_{+y}$ | $J_{-y},H_{-y}$ | $J_{+y},H_{-y}$ | $J_{-y},H_{+y}$ |
|---|---|---|---|---|
| a, $30^\circ$ | a, $30^\circ$ | f, $330^\circ$ | d, $210^\circ$ | c, $150^\circ$ |
| b, $90^\circ$ | a, $30^\circ$ | f, $330^\circ$ | d, $210^\circ$ | c, $150^\circ$ |
| c, $150^\circ$ | a, $30^\circ$ | f, $330^\circ$ | d, $210^\circ$ | c, $150^\circ$ |
| d, $210^\circ$ | a, $30^\circ$ | f, $330^\circ$ | d, $210^\circ$ | c, $150^\circ$ |
| e, $270^\circ$ | a, $30^\circ$ | f, $330^\circ$ | d, $210^\circ$ | c, $150^\circ$ |
| f, $330^\circ$ | a, $30^\circ$ | f, $330^\circ$ | d, $210^\circ$ | c, $150^\circ$ |

These are notebook-recorded states, without a pulse duration, current amplitude,
or numerical trace on this page.

## Page 16 — Six-state sketch

Source: [PDF page 16](../Mn3Sn.pdf#page=16) · [page image](Mn3Sn-pages/page-16.jpg)

### `P16-SIX-STATE-SKETCH`

The relation $J_s=\theta\,\sigma\times J_c$ is annotated $(z,x,y)$.
Six three-spin drawings appear in two rows. $P$ is drawn out of the page and
$I\parallel H$ rightward. No state letters, angles, transition arrows,
or further equations are written on this page. The drawings remain in the
source image without assigning new labels to them.

## Page 17 — Parameters, units, and field response

Source: [PDF page 17](../Mn3Sn.pdf#page=17) · [page image](Mn3Sn-pages/page-17.jpg)

### `P17-FIELD-RESPONSE`

Notes at the top read “$H_x$ small oscillation” and “$H_x$ medium.” Two
qualitative $M_{\mathrm{oct}}$-versus-$t$ sketches show descending traces,
both labeled $J=+x$, with $H_x=0.01$ T and $H_x=-0.01$ T respectively.
The axes have no numerical ticks.

For “$H_x$ large,” the note says:

> it evolves to a new state (m is sensitive to in-plane $H_{\mathrm{eff}}$)

### `P17-PARAMETERS-UNITS`

The constants are written as

$$
\mu_B=9.274\times10^{-24}\ \mathrm{J/T},\qquad
e=1.602\times10^{-19}\ \mathrm C,\qquad
\mathrm T=\mathrm{V\,s/m^2}.
$$

Under source `[3]`, the page writes $M_s=3\mu_B$ and a torque amplitude

$$
T_{ia}=\frac{\hbar J_0\theta_{\mathrm{SH}}}{2e m_s d},
\qquad m_s=6M_s\ \text{per unit cell volume}.
$$

Uppercase and lowercase $M_s/m_s$ are retained as written; their usage changes
within the page. The visible parameter list is:

| Symbol | Notebook value |
|---|---|
| $J$ | 23 meV |
| $D$ | 1.6 meV |
| $k$ | 0.17 meV |
| $\alpha$ | 0.003 |
| $\theta_{\mathrm{SH}}$ | 0.1 |
| $d$ | 40 nm |
| $J_0$ | $6\times10^{11}$ A/m² |
| $H_{\mathrm{ext}}$ | 0.1 T |

These historical notebook values are not the current repository defaults.
Here $J_0$ labels a current density, whereas page 6 uses $J_0$ in an exchange
energy; the notation is not unified in this transcription.

### `P17-EXCHANGE-AND-TORQUE-UNITS`

The exchange and derivative lines read

$$
H=J\sum\vec m_{ia}\cdot\vec m_{ib},\qquad
H_{\mathrm{eff}}=-\frac1{M_s}\frac{\partial H}{\partial\vec m_{ia}}.
$$

The worked exchange-field scale is

$$
H_{\mathrm{ex}}=
\frac{23\times10^{-3}\times1.602\times10^{-19}\ \mathrm J}
{3\times9.274\times10^{-24}\ \mathrm{J/T}}
=132.43\ \mathrm T.
$$

The SOT unit calculation writes

$$
\frac{\hbar J_0\theta_{\mathrm{SH}}}{2eM_s d}
\ \longrightarrow\
\frac{\mathrm{eV\,s}\;\mathrm{A/m^2}}
{e\;\mathrm{A/m}\;\mathrm m}
=\mathrm{V\,s/m^2}=\mathrm T.
$$

The last moment/magnetization lines say

$$
m_s=18\mu_B\quad(\mathrm{J/T}),\qquad
M_s=\frac{m_s}{V}=\frac{\mathrm{J/T}}{\mathrm{m^3}}
=\mathrm{A/m}.
$$

Under source `[5]`, the final lines give

$$
H_{\mathrm{SOT}}=\frac{\gamma\hbar\theta_{\mathrm{SH}}J}{2eM_s d},
\qquad d=10\ \mathrm{nm},\qquad M_s=1.26\ \mathrm T.
$$

Status: `UNRESOLVED`. The notation, gamma factor, and magnetization units in
these different lines are preserved without converting them into a common
field convention or adopting them as code inputs.

## Page 18 — Configuration I under an in-plane field

Source: [PDF page 18](../Mn3Sn.pdf#page=18) · [page image](Mn3Sn-pages/page-18.jpg)

### `P18-IN-PLANE-FIELD`

The heading is “Configuration I,” source `[3]`. Orange spin arrows at green
sites form three-spin motifs. The upper left and upper right drawings have
$H_y$ pointing right, while the lower left drawing has $H_y$ pointing left.
The upper right includes two nearby motifs, one drawn with thicker arrows.
No numerical angles, field amplitudes, time axis, or explicit transition
arrow is provided. The relative arrow directions and thicknesses remain
visible in the source image; no quantitative response is inferred.

## Page 19 — Atomistic model: 5 × 6 example

Source: [PDF page 19](../Mn3Sn.pdf#page=19) · [page image](Mn3Sn-pages/page-19.jpg)

### `P19-ATOMISTIC-5X6`

This landscape page contains an annotated photograph headed “Version 3,”
“Atomistic model — zhengde,” “5 × 6 example,” and “最近邻” (“nearest neighbors”).
Black hexagons are overlaid with blue circled sites, red crossed circles,
green/orange tracing lines, and row/column indices. $L$ points right;
$W$ points up, while row labels increase downward.

Red note, connected to the opposing W/row directions:

> their opposite may produce error

Rows 1–5 and column numbers are visible. An orange-circled crossed site at
the right is labeled 7. A green “repeat” arrow points downward.

### `P19-NEIGHBOR-DIRECTION-SKETCHES`

The right-hand notes use four row classes written as division-by-four
remainders “1,” “2,” “3,” and “4,” each split into 奇 / 偶 (odd / even columns).
The sketch directions, read in the image frame, are:

| Written row class | Odd column sketch | Even column sketch |
|---|---|---|
| 1 | up, up-right, down, down-right | 0 |
| 2 | up, left, right, down-left | up-left, left, right, down |
| 3 | 0 | up, up-right, down, down-right |
| 4 | up-left, left, right, down | up, left, right, down-left |

The bottom shorthand reads approximately “行 /4 … 1 … 3? 奇” and
“列 … 奇”; the intervening marks are `[unclear]`. The four-class sketches
are not converted into array offsets, because the page explicitly questions
the relation between physical W and row direction.

## Page 20 — Lattice and atom-type array

Source: [PDF page 20](../Mn3Sn.pdf#page=20) · [page image](Mn3Sn-pages/page-20.jpg)

### `P20-OCCUPATION`

The top lattice drawing labels rows 1–5 and columns 1–6, with green circled
sites and orange spin arrows. A larger lattice drawing underneath repeats
the connected triangles without spin arrows.

The bottom “Atom-type” table is transcribed exactly in the displayed order:

```text
        c1 c2 c3 c4 c5 c6
row 1    1  3  1  3  1  3
row 2    1  1  1  1  1  1
row 3    3  1  3  1  3  3
row 4    1  1  1  1  1  1
row 5    1  3  1  3  1  3
```

In particular, row 3 column 6 is written **3**, not a continuation of the
alternating 3/1 pattern. That boundary entry is retained. The notebook image
does not supply a separate numeric legend on this page; the repository's
documented 1/3 occupancy convention is described separately in the README.

## Page 21 — Anisotropy and exchange masks

Source: [PDF page 21](../Mn3Sn.pdf#page=21) · [page image](Mn3Sn-pages/page-21.jpg)

### `P21-ANISOTROPY-EXCHANGE-MASKS`

The heading is “zhengde paper / Configuration II.” The triad has $z$ up,
$y$ right, and $x$ out of the page. The spin sketch is f, $330^\circ$;
$P$ is out of the page and $I\parallel H$ is rightward. An orange axis
sketch marks $\pi/6$ and $5\pi/6$.

Each array below has five rows and six columns. `X` preserves red/blue crosses,
`A` the orange coupling symbol, and `0` the orange zero. The page does not
give a numerical value or unit for `A`. The anisotropy table is labeled
“anisotropy — y, z”:

```text
3π/2   X    3π/2   X    3π/2   X
5π/6  π/6   5π/6  π/6   5π/6  π/6
 X    3π/2   X    3π/2   X     X
π/6   5π/6  π/6   5π/6  π/6   5π/6
3π/2   X    3π/2   X    3π/2   X
```

### `P21-L-MASKS`

`A_next L`:

```text
0 X 0 X 0 X
A A A A A 0
X 0 X 0 X X
A A A A A 0
0 X 0 X 0 X
```

`A_prev L`:

```text
0 X 0 X 0 X
0 A A A A A
X 0 X 0 X X
0 A A A A A
0 X 0 X 0 X
```

Margin notes say odd rows are 0 and even rows are 1; the next-L last column
and previous-L first column are 0. The matrices themselves use `A`, not `1`.
Both notations are retained without assigning a numerical normalization.

### `P21-W-MASKS`

`A_next W`:

```text
0 X 0 X 0 X
A 0 A 0 A 0
X A X A X X
0 A 0 A 0 0
A X A X A X
```

`A_prev W`:

```text
A X A X A X
0 A 0 A 0 0
X A X A X X
A 0 A 0 A 0
0 X 0 X 0 X
```

Next-W notes: first row all 0; second-row odd columns 1; odd rows all 1;
fourth-row even columns 1. Previous-W notes: odd rows all A; second-row
even columns 1; fourth-row odd columns 1; last row all 0. Crossed slots and
the displayed edge exceptions remain as drawn.

### `P21-DIAGONAL-LEFT-MASKS`

`A_next W left`:

```text
0 X 0 X 0 X
0 A 0 A 0 A
X 0 X 0 X X
0 0 A 0 A 0
0 X 0 X 0 X
```

`A_prev W left`:

```text
0 X 0 X 0 X
0 0 A 0 A 0
X 0 X 0 X X
0 A 0 A 0 A
0 X 0 X 0 X
```

Next-W-left notes: first column and first row 0; row 2 even columns A;
odd rows 0; row 4 odd columns A. Previous-W-left notes: first column 0;
row 2 odd columns have “prev left”; odd rows have no “prev left”;
row 4 even columns have “prev left”; last row 0.

### `P21-DIAGONAL-RIGHT-MASKS`

`A_next W right`:

```text
0 X 0 X 0 X
0 0 0 0 0 0
X A X A X X
0 0 0 0 0 0
A X A X A X
```

`A_prev W right`:

```text
A X A X A X
0 0 0 0 0 0
X A X A X X
0 0 0 0 0 0
0 X 0 X 0 X
```

Next-W-right notes: last column 0, first row 0, odd rows 1, even rows 0.
Previous-W-right notes: last column 0, odd rows A, even rows 0,
row 3 A, last row 0. The written “last column 0” notes coexist with crossed
slots in that column; the crosses are preserved.

Status: these are literal handwritten arrays and margin rules. Their edge
behavior, reciprocity, and agreement with MATLAB have not been checked in this
extraction. Page 23 repeats the arrays with a different selection of annotations.

## Page 22 — Atomistic model: 7 × 6 example

Source: [PDF page 22](../Mn3Sn.pdf#page=22) · [page image](Mn3Sn-pages/page-22.jpg)

### `P22-ATOMISTIC-7X6`

This landscape page repeats the Version 3 nearest-neighbor photograph from
page 19 but labels it “7 × 6 example.” Blue site circles and green/orange
tracing lines extend through two additional rows. Right-edge sites near the
third and seventh rows are circled and crossed in orange.

$L$ still points right, $W$ up, and the row labels increase down. The red note
“their opposite may produce error,” the green “repeat” arrow, and the
four row-class odd/even neighbor-direction sketches remain. The directional
table is the same as on page 19:

| Written row class | Odd column sketch | Even column sketch |
|---|---|---|
| 1 | up, up-right, down, down-right | 0 |
| 2 | up, left, right, down-left | up-left, left, right, down |
| 3 | 0 | up, up-right, down, down-right |
| 4 | up-left, left, right, down | up, left, right, down-left |

The lower row/column shorthand is again `[unclear]`. This page is retained
separately because it extends the example; it is not discarded as a duplicate.

## Page 23 — References

Source: [PDF page 23](../Mn3Sn.pdf#page=23) · [page image](Mn3Sn-pages/page-24.jpg)

### `P24-REFERENCES`

The heading is “Refs.” The five abbreviated entries read:

| Notebook label | Year | Journal as written | Name as written |
|---|---|---|---|
| [1] | 2017 | PRB | Suzuki |
| [2] | 2015 | nature | Satoru Nakatsuji |
| [3] | 2020 | nature | Hanshen Tsai |
| [4] | 2022 | nature | Tomoya Higo |
| [5] | 2021 | nat mat | Yutaro Takeuchi |

These are literal notebook labels, not independently resolved bibliographic
entries. Article titles, full author lists, volume/article numbers, DOIs, and
exact source locations are **REFERENCE MISSING**. The spelling of “Hanshen”
is retained from the handwritten entry. “Zhengde paper” elsewhere in the
notebook is not explicitly assigned one of these five numbers.
