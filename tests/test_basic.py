def test_data_loader():
    from holoseq.data import load
    path = "data/mUroPar1H1H2.paf_cisH1_hseq.gz"
    assert load(path)[0] is not None