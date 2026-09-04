import math

from dockstudio.core import utils


def test_crc32_seed_reproducible():
    a = utils.crc32_seed("rec1", "lig1")
    b = utils.crc32_seed("rec1", "lig1")
    assert a == b
    assert 0 <= a <= 0x7FFFFFFF
    c = utils.crc32_seed("rec1", "lig2")
    assert a != c


def test_safe_name():
    assert utils.safe_name("a/b c:分子") == "a_b_c"
    assert utils.safe_name("4DFR") == "4DFR"


def test_ascii_only():
    assert utils.ascii_only("中文路径") == ""


def test_distance():
    assert utils.distance([0, 0, 0], [3, 4, 0]) == 5.0
