from pdb2mcsm.core import parse_pdb, assign_chain_ids, write_clean, build_ala_scan

def test_two_chains_via_ter():
    header, chains, warnings, method = parse_pdb("tests/fixtures/mini_raw.pdb")
    assign_chain_ids(chains)
    assert method == "TER records"          # TER로 나뉘었어야
    assert len(chains) == 2                  # 두 사슬
    assert warnings == []                    # 문제 없음
    assert chains[0].residues[2][1] == "HIS" # HIE가 HIS로 정규화됐는지 (첫 사슬 3번째 잔기)

def test_clean_output_has_no_water_or_h(tmp_path):
    header, chains, warnings, method = parse_pdb("tests/fixtures/mini_raw.pdb")
    assign_chain_ids(chains)
    out = write_clean(str(tmp_path/"clean.pdb"), header, chains)
    text = open(out).read()
    assert "HOH" not in text and " H1 " not in text   # 물·수소 제거됨
    assert "HIS" in text and "HIE" not in text         # 정규화됨

def test_ala_scan_format():
    header, chains, warnings, method = parse_pdb("tests/fixtures/mini_raw.pdb")
    assign_chain_ids(chains)
    muts = build_ala_scan(chains, "A", 1, 3)
    assert muts == ["A S1A", "A R2A", "A H3A"]