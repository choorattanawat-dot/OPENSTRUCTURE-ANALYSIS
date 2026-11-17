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
2. P-Delta, stability coefficient Q, story drift, ELF/RSA ตาม มยผ. 1301/1302-61
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
