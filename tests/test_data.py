from pathlib import Path

from holoseq.data import load

TEST_DIR = Path(__file__).parent.resolve()
DATA_DIR = TEST_DIR.joinpath("../data")


def test_load():
    path = DATA_DIR / "mUroPar1H1H2.paf_cisH1_hseq.gz"
    assert load(path)[0] is not None
