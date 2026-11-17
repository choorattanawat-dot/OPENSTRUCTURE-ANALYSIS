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
| T1 (s) | 2.697 | 2.718 | +0.8 % |
| T2 (s) | 1.784 | 1.786 | +0.1 % |
| T3 (s) | 1.455 | 1.462 | +0.5 % |
| T4 (s) | 0.732 | 0.730 | −0.3 % |

Mode 5–12 ต่างไม่เกิน 3.1 % และ participating mass ratio ของทุก mode ตรงกันภายใน 0.02

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

## ข้อตกลงการแปลง (ห้ามเปลี่ยนเงียบ ๆ ถ้าเปลี่ยนต้องขึ้นเวอร์ชันและรัน test ใหม่)

- สร้างเฉพาะ object ที่มี LINE/AREA ASSIGNS
- Mesh: กริดโครงสร้างเดียวผ่านทุกพิกัดสำคัญ แบ่งไม่เกิน 1.25 m (พื้น) และ 0.5 m (ผนัง) ต้องเป็นขอบขนานแกน
- พื้น/ผนัง: `ShellMITC4` + `ElasticMembranePlateSection` ใช้ modifier F11 (membrane) และ M11 (bending)
- **ถ้า M12MOD ของพื้นต่างจาก M11MOD** (เช่น ตั้ง M11 = M22 = 0.25 แต่ปล่อย M12 = 1 ตามค่าเริ่มต้น ETABS)
  จะใช้ `LayeredShell` + `ElasticOrthotropic` เพื่อให้ twisting stiffness ตรงกับ ETABS
  (กรณีอาคารอ้างอิง: ถ้าไม่ทำ T1 จะเพี้ยน +14 %)
- เสา/คาน: `elasticBeamColumn` (Euler-Bernoulli), modifier จาก section x assign, แบ่งที่ทุก node บนชิ้นส่วน
- พื้นซ้อนทับกัน: ชิ้นแรกให้ความแข็ง/น้ำหนัก ส่วน area load ของทุกชิ้นที่ทับถูกใส่ครบ (เหมือน ETABS)
- มวล: self-mass ของ element เท่านั้น เฉพาะ UX/UY รวมที่ระดับชั้น
- ฐานรองรับ: จุดที่มี RESTRAINT ที่ Base และ node ฐานผนังทุกตัว (UX UY UZ)

**ยังไม่ได้จำลอง:** rigid end zone, cardinal point offset, panel zone, link, tendon/PT balancing load,
แรงลม, shear deformation ของ frame, หน้าตัดเหล็ก, ผนัง/พื้นที่ไม่ขนานแกน

## แผนงาน

1. ✅ ตัวสร้างโมเดล + ตรวจสอบกับ ETABS (gravity, mass, modal)
2. ✅ P-Delta, stability Q, story drift ⬜ ELF/RSA ตาม มยผ. 1301/1302-61
3. ตรวจออกแบบ คสล. (ACI 318 / วสท.)
4. Nonlinear: fiber section, pushover, NLTHA, IDA
5. รายงาน + 3D viewer

## ทดสอบ

```
pytest                                   # unit test (ไม่ต้องใช้โมเดลโครงการ)
OPENSTRUCTURE_REF_E2K=... OPENSTRUCTURE_REF_TABLES=... pytest    # เทียบ ETABS กับโมเดลจริง
```

ไฟล์โมเดลโครงการ (`.e2k`) ไม่เก็บใน repo นี้

OpenSees ใช้ฟรีสำหรับงานวิจัย การศึกษา และใช้ภายใน การขายซอฟต์แวร์ที่มี OpenSees ต้องขอ licence จาก UC Berkeley
