"""DPT 1301/1302-61 values checked against the standard (Khon Kaen, Mueang district)."""
import pytest

from openstructure.dpt import Site, k_exponent, period_approx


def test_khon_kaen_site_d():
    s = Site.named("ขอนแก่น อ.เมือง", soil="D")
    assert (s.Fa, s.Fv) == (1.6, 2.4)
    assert s.SDS == pytest.approx(0.05653, abs=1e-5)
    assert s.SD1 == pytest.approx(0.0480, abs=1e-5)
    assert s.Ts == pytest.approx(0.849, abs=1e-3)
    assert s.sdc() == "ก"
    assert s.sa_elf(0.5) == pytest.approx(s.SDS)
    assert s.sa_elf(2.0) == pytest.approx(s.SD1 / 2)
    assert s.sa_dynamic(0.0) == pytest.approx(0.4 * s.SDS)


def test_khon_kaen_site_e_exception():
    s = Site.named("ขอนแก่น อ.เมือง", soil="E")
    assert s.SDS == pytest.approx(0.0883, abs=1e-4)
    assert s.SD1 == pytest.approx(0.070, abs=1e-4)
    assert s.sdc() == "ข"              # Table 1.6-2 governs
    assert s.sdc(Ta=0.5) == "ก"        # Ta < 0.8 Ts -> Table 1.6-1 alone


def test_period_and_k():
    assert period_approx(30) == pytest.approx(0.6)
    assert k_exponent(0.4) == 1 and k_exponent(3) == 2 and k_exponent(1.5) == pytest.approx(1.5)
