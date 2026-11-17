# OPENSTRUCTURE-ANALYSIS

เครื่องมือวิเคราะห์โครงสร้างของเราเอง บน **OpenSeesPy** อ่านโมเดลจากไฟล์ข้อความ ETABS (`.e2k`)
แล้วสร้างโมเดล 3 มิติ วิเคราะห์ และเทียบผลกับ ETABS โดยไม่ต้องใช้ ETABS API

หน่วยภายในโมเดล OpenSees: **kN, m, s, t** ไฟล์ `.e2k` เป็นหน่วยใดก็ได้ (kgf-m, kN-m, N-mm, tonf-m, lb-in, kip-ft)

## ติดตั้ง

```
pip install -e .          # ต้องมี Python 3.10+ (Windows: openseespy ต้องใช้ Python 3.12)
```

## ใช้งาน

```python
from openstructure import Model

m = Model.from_e2k("building.e2k")
m.info()            # จำนวน node / shell / frame, คำเตือน, สิ่งที่ยังไม่ได้จำลอง
m.total_mass        # มวลรวม (t)
m.gravity()         # base reaction ของแต่ละ load pattern
m.modal(12)         # คาบธรรมชาติ + modal participating mass ratio
m.lateral_share()   # สัดส่วนแรงเฉือน ผนัง/เสา, story drift, torsion ratio
m.story_drift("X")  # story drift, drift ratio, torsion ratio (แรงด้านข้างกำหนดเองได้)
m.stability()       # ค่า Q ต่อชั้น + P-Delta amplification
m.write_tcl("model.tcl")   # คำสั่ง OpenSees ทั้งหมด ไว้ตรวจสอบย้อนหลัง
```

เทียบกับ ETABS (export ตารางผลจาก ETABS เป็น text):

```
python -m openstructure building.e2k etabs_tables.e2k
```

ตารางที่ใช้: `Base Reactions`, `Modal Participating Mass Ratios`, `Program Control`

## ผลตรวจสอบกับ ETABS (อาคารอ้างอิง, ETABS 23.1.1)

อาคาร คสล. เสา 40x80, พื้น PT flat slab 250 mm, core wall, 9 ชั้น, ไฟล์หน่วย kgf-m

| รายการ | ETABS | OpenSees | ต่าง |
|---|---|---|---|
| Base FZ DL (kN) | 34480.6 | 34489.6 | +0.03 % |
| Base FZ SDL (kN) | 1874.2 | 1874.2 | 0.00 % |
| Base FZ LL (kN) | 8982.0 | 8982.0 | 0.00 % |
| T1 (s) | 2.697 | 2.721 | +0.9 % |
| T2 (s) | 1.784 | 1.789 | +0.3 % |
| T3 (s) | 1.455 | 1.463 | +0.5 % |
| T4 (s) | 0.732 | 0.731 | −0.1 % |

Mode 5–12 ต่างไม่เกิน 4 % และ participating mass ratio ของทุก mode ตรงกันภายใน 0.02

แรงในชิ้นส่วน: แรงอัดเสาและผนัง core ตรงกับ ETABS ภายใน 1–2 % ส่วนแรงลมในผนัง core ได้แรงเฉือนน้อยกว่า ETABS ราว 5 %
และโมเมนต์น้อยกว่าราว 10 % (การแบ่งแรงระหว่างผนังกับโครงพื้น-เสาไวต่อการจำลองจุดต่อพื้น-เสา ยังหาสาเหตุไม่จบ)

## Story drift, ค่า Q และ P-Delta

```python
from openstructure.stability import format_stability
print(format_stability(m.stability()))                       # 1.2D + 1.0L, แรง F ~ m*z รวม 1000 kN
m.stability(gravity_factors={"DL": 1, "SDL": 1, "LL": 1},    # น้ำหนักไม่คูณ factor = theta ของ มยผ.
            forces={"X": {"FL.7": 50.6, ...}, "Y": {...}})   # ใส่แรงลม/แรงแผ่นดินไหวต่อชั้นเองได้ (kN)
```

หรือ `python -m openstructure building.e2k --stability`

- Story drift = ค่าเฉลี่ยถ่วงมวลของ drift ทุกแนวดิ่งที่มี node ทั้งสองระดับ, torsion ratio = d_max / d_avg
  (มยผ. 1301/1302-61 ข้อ 2.4.2.1: > 1.2 แบบ 1ก, > 1.4 แบบ 1ข)
- Q = ΣP·Δ / (V·h) ตาม ACI 318-14 ข้อ 6.6.4.3 (Q > 0.05 = sway) ถ้าใช้น้ำหนักไม่คูณ factor จะเท่ากับ
  θ ของ มยผ. 1301/1302-61 สมการ 3.8-1 (θ ≤ 0.10 ไม่ต้องคิด P-Delta, θmax = 0.5/(β·Cd) ≤ 0.25)
- P-Delta: เสาใช้ `geomTransf PDelta` ภายใต้ gravity คงที่ แล้วรายงาน drift ที่ขยายเทียบกับ 1/(1−Q)
  ผนังและพื้นเป็น shell ไม่มี geometric stiffness จึงได้ P-Delta เฉพาะแรงอัดที่ลงเสา
- ตรวจกับสูตรมือ (เสายื่น): drift = FL³/3EI และ amplification = 1/(1 − PL²/3EI) ตรงกัน

## แรงแผ่นดินไหว มยผ. 1301/1302-61 (นอกแอ่งกรุงเทพ)

```python
from openstructure.dpt import Site, elf, format_elf
site = Site.named("ขอนแก่น อ.เมือง", soil="D", occupancy="II")   # หรือ Site(Ss, S1, soil, occupancy)
site.SDS, site.SD1, site.sdc()                                    # 0.0565, 0.0480, "ก"
T = m.modal(3)["periods_s"]
res = elf(m, site, system="bf_ordinary_rc_wall", T_model={"X": T[2], "Y": T[1]})
print(format_elf(res))
```

ได้ Fa/Fv (ตาราง 1.4-2, 1.4-3), SDS/SD1, สเปกตรัม (รูป 1.4-1 ถึง 1.4-4), ประเภทการออกแบบ (ตาราง 1.6-1, 1.6-2),
Ta = 0.02H, Cs = Sa·I/R ≥ 0.01, การกระจายแรง Cvx, δx = Cd·δxe/I เทียบ Δa (ตาราง 2.11-1) และ θ (สมการ 3.8-1)
ยังไม่มี: สเปกตรัมแอ่งกรุงเทพ, accidental torsion, RSA

## ตรวจออกแบบเสา คสล. (ACI 318-14)

```python
from openstructure.design.column import check_columns, format_columns
rows, secs = check_columns(m, lu_deduct=0.25)          # ใช้ load combination และเหล็กเสริมจาก e2k
rows, secs = check_columns(m, fc=37.3, combos={"U2": {"DL": 1.2, "SDL": 1.2, "LL": 1.6}})
print(format_columns(rows, secs, top=20))               # เรียงจาก D/C มากไปน้อย
m.frame_forces()                                         # แรงในเสา/คานของแต่ละ load pattern
```

- แรง Pu, M2, M3 มาจาก OpenSees (รวม pattern ตาม combo แบบเชิงเส้น) แรงอัดเสาตรงกับ ETABS ภายใน 1–2 %
- กำลังหน้าตัด: interaction surface P-M2-M3 จาก strain compatibility (εcu 0.003, Whitney block, φ 0.65–0.90
  ตาม εt, จำกัด φPn ≤ 0.80φPo) แล้วหา D/C ตามแนวรังสีจากจุดกำเนิด
- ความชะลูด: k = 1, r = 0.3h, (EI)eff = 0.4EcIg/(1+βdns), Cm, M2,min = Pu(15+0.03h), δns = Cm/(1−Pu/0.75Pc)
  ชั้น sway ใช้ขีดจำกัด klu/r ≤ 22 (ค่าเริ่มต้น `sway=True`) ชั้น nonsway ใช้ 34+12(M1/M2) ≤ 40
- `PATTERN "R-a-b"` อ่านเป็น a เส้นต่อด้านที่ขนานแกน 3 และ b เส้นต่อด้านที่ขนานแกน 2 ถ้าแบบจริงต่างให้ใส่
  `n_face3=`, `n_face2=`
- ยังไม่ได้ทำ: δs สำหรับ combo แรงด้านข้าง (ให้ใช้ผล P-Delta), live load reduction, แรงเฉือน, รายละเอียดเหล็กปลอก

## ตรวจออกแบบ lift core / ผนัง (ACI 318-14)

```python
from openstructure.etabs_tables import story_forces
from openstructure.design.wall import check_walls, format_walls
W = story_forces("etabs_tables.e2k", "W", bottom_story="FL.2")          # แรงลมต่อชั้นจาก ETABS
m.add_lateral_case("WX", "X", W[1]["X"], W[1]["y_of_X"])                # ใส่เป็น load case ของเรา
m.add_lateral_case("WY", "Y", W[2]["Y"], W[2]["x_of_Y"])
combos = {"0.9D+1.6WY": {"DL": 0.9, "SDL": 0.9, "WY": 1.6}, ...}
rows, secs = check_walls(m, combos, fc=23.5, rho_v=0.0025, rho_h=0.0025)
print(format_walls(rows, secs))
m.wall_forces()      # P, Vx, Vy, Mx, My, T ของ pier ที่โคนผนังแต่ละชั้น
```

- ตัดหน้าตัดแนวราบผ่านผนังทุกแผ่นของ pier ที่โคนชั้น (รูปกล่อง, C, L ได้หมด) แล้วรวมแรงรอบ centroid
- P-M-M: strain compatibility ของหน้าตัดจริง เหล็กยืนกระจายตาม ρv, φ 0.65–0.90, φPn ≤ 0.80φPo
- แรงเฉือน (ข้อ 11.5.4): Vn = (0.17√f'c + ρt·fy)·Σ(t·0.8lw) ≤ 0.83√f'c·Σ(t·0.8lw), φ 0.75
- e2k ไม่มีเหล็กผนัง ต้องใส่ ρv, ρh เอง (ค่าเริ่มต้น 0.0025)
- ยังไม่ได้ตรวจ: boundary element, ความชะลูดผนัง, คานทับหลัง (lintel), แรงบิด

## Base reaction รายจุด

```python
from openstructure.reactions import base_reactions, format_reactions
res = base_reactions(m, combos=m.M["combos"])        # ทุก load pattern + combo
print(format_reactions(res))                          # ตาราง Fz ต่อจุด (เลือก component="Fx" ฯลฯ ได้)
```

รายงานตาม label จุดของ ETABS ส่วน node ฐานผนังที่ไม่ใช่จุดของ ETABS รวมเข้าจุดที่ใกล้ที่สุด และบวกน้ำหนักที่ลงบนจุดรองรับ
โดยตรงกลับเข้า Fz ให้ครบ เทียบ ETABS: เสาตรงภายใน 2 % (ต้นมุม 4 %), ผลรวม core ตรงภายใน 1 %

## Pile cap (ACI 318-14 ข้อ 13.4)

```python
from openstructure.design.pilecap import PileCap, pile_layout, design_pilecap, format_pilecap
cap = PileCap(pile_layout(10, 1.05, rotate=True), pile_size=0.35, h=1.20, col=(0.4, 0.8), fc=23.5, fy=390)
r = cap.check(Pu=7042)                 # flexure, one-way shear, punching ที่เสาและที่เข็ม
cap, res = design_pilecap(10, 0.35, (0.4, 0.8), (4813, 0, 0), [("U1", 7042, 0, 0)], fc=23.5)   # หาความหนาน้อยสุด
print(format_pilecap(cap, res))
```

- cap แข็ง เข็มรับแรงแกนอย่างเดียว: R = P/n + Mx·y/Σy² − My·x/Σx²; ผังเข็มมาตรฐาน 1–12 ต้น (10 ต้น = 3-4-3)
- แรงเฉือนนับแรงเข็มบางส่วนตามข้อ 13.4.2.5; ยังไม่มี strut-and-tie (13.2.6.3), แรงด้านข้างเข็ม, เข็มรับแรงถอน

## ข้อตกลงการแปลง (ห้ามเปลี่ยนเงียบ ๆ ถ้าเปลี่ยนต้องขึ้นเวอร์ชันและรัน test ใหม่)

- สร้างเฉพาะ object ที่มี LINE/AREA ASSIGNS
- Mesh: กริดโครงสร้างเดียวผ่านทุกพิกัดสำคัญ แบ่งไม่เกิน 1.25 m (พื้น) และ 0.5 m (ผนัง) ต้องเป็นขอบขนานแกน
- พื้น/ผนัง: `ShellMITC4` + `ElasticMembranePlateSection` ใช้ modifier F11 (membrane) และ M11 (bending)
- **ถ้า modifier แรงเฉือนของพื้น (F12, M12) ต่างจากค่าแนวตรง (F11, M11)** เช่น ตั้ง M11 = M22 = 0.25 แต่ปล่อย
  M12 = 1 ตามค่าเริ่มต้น ETABS จะใช้ `LayeredShell` + `ElasticOrthotropic` แบบแกนกลาง + ผิวสองด้าน แล้วแก้สมการหา
  E และ G12 ของแต่ละชั้นให้ผลรวมตรงกับ F11, M11, F12, M12 พอดี (กรณีอาคารอ้างอิง: ถ้าไม่ทำ T1 จะเพี้ยน +14 %)
  ส่วน transverse shear ไม่คูณ modifier เพราะ ShellThin ของ ETABS ไม่มี shear deformation
- เสา/คาน: `elasticBeamColumn` (Euler-Bernoulli), modifier จาก section x assign, แบ่งที่ทุก node บนชิ้นส่วน
- พื้นซ้อนทับกัน: ชิ้นแรกให้ความแข็ง/น้ำหนัก ส่วน area load ของทุกชิ้นที่ทับถูกใส่ครบ (เหมือน ETABS)
- มวล: self-mass ของ element เท่านั้น เฉพาะ UX/UY รวมที่ระดับชั้น
- ฐานรองรับ: จุดที่มี RESTRAINT ที่ Base และ node ฐานผนังทุกตัว (UX UY UZ)

**ยังไม่ได้จำลอง:** rigid end zone, cardinal point offset, panel zone, link, tendon/PT balancing load,
แรงลม, shear deformation ของ frame, หน้าตัดเหล็ก, ผนัง/พื้นที่ไม่ขนานแกน

## แผนงาน

1. ✅ ตัวสร้างโมเดล + ตรวจสอบกับ ETABS (gravity, mass, modal)
2. ✅ P-Delta, stability Q, story drift, ELF ตาม มยผ. 1301/1302-61 ⬜ RSA, accidental torsion, แอ่งกรุงเทพ
3. ตรวจออกแบบ คสล. (ACI 318-14): ✅ เสา ✅ lift core ⬜ พื้น PT, บันได, ฐานราก
4. Nonlinear: fiber section, pushover, NLTHA, IDA
5. รายงาน + 3D viewer

## ทดสอบ

```
pytest                                   # unit test (ไม่ต้องใช้โมเดลโครงการ)
OPENSTRUCTURE_REF_E2K=... OPENSTRUCTURE_REF_TABLES=... pytest    # เทียบ ETABS กับโมเดลจริง
```

ไฟล์โมเดลโครงการ (`.e2k`) ไม่เก็บใน repo นี้

OpenSees ใช้ฟรีสำหรับงานวิจัย การศึกษา และใช้ภายใน การขายซอฟต์แวร์ที่มี OpenSees ต้องขอ licence จาก UC Berkeley
