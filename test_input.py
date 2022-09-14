from openfold import format_input
import pytest

def test_format_input_1():
    assert format_input("MaAHKGAEHHHKAAEHHEQAAKHHHAAAEHHeKGEHEQAAHHADTAYAHHKHAeeHAAQAAKHDAEHHAPKPH") == "MAAHKGAEHHHKAAEHHEQAAKHHHAAAEHHEKGEHEQAAHHADTAYAHHKHAEEHAAQAAKHDAEHHAPKPH"

def test_format_input_2():
    with pytest.raises(Exception):
        format_input("xyz")