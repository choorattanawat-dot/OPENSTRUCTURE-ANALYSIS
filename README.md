# OPENSTRUCTURE-ANALYSIS

Structural analysis and reinforced-concrete member design for 2D frames. Version 0.1.

## สรุปภาษาไทย

ไลบรารีภาษาไพธอนสำหรับวิเคราะห์โครงข้อแข็ง 2 มิติแบบสถิตเชิงเส้น ด้วย OpenSeesPy และออกแบบหน้าตัดคอนกรีตเสริมเหล็กโดยวิธีกำลัง ตามมาตรฐาน วสท. 1008 และ ACI 318-19 ใช้ได้จากบรรทัดคำสั่งและจาก Jupyter

แกนคำนวณแยกจากส่วนอ่านไฟล์และส่วนแสดงผล เพื่อให้ต่อไปเพิ่มการวิเคราะห์ไม่เชิงเส้น เช่น pushover หรือไฟเบอร์ ได้โดยไม่ต้องเขียนโมเดลใหม่ หน่วยเป็นระบบ SI และระบุไว้ในผลลัพธ์ รายงานการออกแบบแสดงทีละขั้นเพื่อให้วิศวกรตรวจทานได้

ข้อสมมุติของตัวคูณและสูตรที่ยังต้องให้วิศวกรยืนยันอยู่ใน [docs/CODE_ASSUMPTIONS.md](docs/CODE_ASSUMPTIONS.md)

## What v0.1 does

- 2D frame model: nodes, supports, rectangular members, nodal loads, member point loads, and uniform loads. Load cases and factored combinations, including ACI 318-19 §5.3.1 and the EIT-era gravity combination 1.4D + 1.7L.
- Linear-static analysis with OpenSeesPy: displacements, reactions, member end forces, and N–V–M diagrams along each member.
- RC design, side by side:
  - beam flexure, singly or doubly reinforced
  - beam shear, vertical stirrups
  - short-column axial plus uniaxial bending
- YAML or JSON input and a command-line runner.
- Tests against classical beam theory and hand-worked design numbers. Sources are in [docs/BENCHMARKS.md](docs/BENCHMARKS.md).

The analysis model does not know which design code will be used later. Nonlinear analysis is not implemented. The model, the load combinations, and the section objects are the parts a later pushover or fiber-section engine should reuse.

## Install

Python 3.10 or newer.

```bash
python -m pip install -e ".[dev]"
```

OpenSeesPy is required for analysis. Design calculations import without calling it, but it is a package dependency.

## Command line

```bash
openstructure run examples/two_bay_two_storey.yaml --out results
```

That writes `results/results.json`, `results/report.md`, and force diagrams under `results/diagrams/`. The example is a 2-bay, 2-storey frame. For the EIT gravity combination 1.4D + 1.7L the first-floor beam B10 has a sagging moment of about 96 kN·m and a hogging moment of about 186 kN·m. Both codes ask for about 1282 mm² of tension steel for that hogging moment. Shear at d beyond the column face is about 142 kN, detailed as 2-leg 9 mm stirrups at 200 mm. The base column, with 8-DB20, is within the P–M diagram (utilization about 0.28). These are example loads, not a project takeoff.

A notebook of the same frame is `examples/two_bay_two_storey.ipynb`.

## Units

| Quantity | Input and reports | Inside the solver |
|---|---|---|
| Coordinates | m | mm |
| Section size | mm | mm |
| Force | kN | N |
| Uniform load | kN/m | N/mm (the same number) |
| Moment | kN·m | N·mm |
| Stress | MPa | MPa = N/mm² |

Global axes: x to the right, y upward, moment counter-clockwise.

Member local axes: x from node i to node j, y at 90° counter-clockwise from x. For a beam drawn left to right, local y is up. For a column drawn upward, local y is horizontal, in the negative global x direction. On a square column the inertia is the same either way. `h_mm` is the section dimension along local y (the bending depth). `b_mm` is the width perpendicular to the frame.

Section forces, and the diagrams:

- Axial force positive in tension.
- Shear positive when the left side of a cut pushes the right side toward local +y. A downward-loaded simple beam has positive shear at the left end.
- Moment positive when it compresses the local +y face. Sagging is positive on a left-to-right beam. Hogging is negative.

Joint end forces stored in the JSON file are the forces the joint exerts on the member, in local axes. The report's `stations` are the internal N, V, and M above.

## Design codes

```python
from openstructure.design import ACI318, EIT1008, design_beam_flexure

aci = design_beam_flexure(200, b_mm=300, d_mm=450, h_mm=500, fc_mpa=28, fy_mpa=400, code=ACI318())
eit = design_beam_flexure(200, b_mm=300, d_mm=450, h_mm=500, fc_mpa=28, fy_mpa=400, code=EIT1008())
print(aci.report.to_markdown())
```

`get_code("ACI318-19")` and `get_code("EIT1008")` build the same objects. A new code, including a later steel code, implements the same methods (`beta1`, `phi_column`, `vc_n`, `as_min_beam_mm2`, and the combination preset) and is selected by name. Nothing in the frame model imports a design code.

Read [docs/CODE_ASSUMPTIONS.md](docs/CODE_ASSUMPTIONS.md) before using a result on a real member. In particular, the EIT wind and earthquake factors, the EIT column φ transition, ρmax = 0.75 ρb, and the ksc shear coefficients are assumptions flagged for confirmation. They are repeated at the bottom of every calculation report.

## Smallest analysis

```python
from openstructure import LinearStaticEngine, Model

model = Model()
model.add_material("C28", ec_mpa=24870, fc_mpa=28, fy_mpa=400)
model.add_rectangular_section("B30x50", "C28", b_mm=300, h_mm=500)
model.add_node("A", 0, 0)
model.add_node("B", 6, 0)
model.add_support("A", pinned=True)
model.add_support("B", roller_x=True)
model.add_member("AB", "A", "B", "B30x50")
model.add_load_case("DL", kind="dead")
model.add_member_udl_global("AB", "DL", wy_kN_per_m=-10)
model.add_combination("DL", {"DL": 1.0})

result = LinearStaticEngine().analyze(model)
beam = result.by_combination["DL"].members["AB"]
print(max(station.m_kNm for station in beam.stations))  # 45 kN·m = w L^2 / 8
```

`AnalysisOptions(kind="linear_static", geometric_transformation="Linear")` is the verified path. `PDelta` can be requested on the same model. Any other `kind` raises an error that names what is not implemented yet, instead of silently running a linear solution.

## Development

```bash
ruff check src tests
pytest
```

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT. OpenSeesPy is a separate package with its own license.
