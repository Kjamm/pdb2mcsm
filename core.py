from dataclasses import dataclass, field

THREE_TO_ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
                "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
                "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

@dataclass
class Atom :
    record : str
    name : str
    resname : str
    chain_id : str
    resnum : int
    icode : str
    element : str
    raw : str

def read_atoms(path) :
    atoms = []

    with open(path) as fh :
        for line in fh :
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM") :
                continue
            atoms.append(Atom(
                record=rec,
                name=line[12:16].strip(),
                resname=line[17:20].strip(),
                chain_id=line[21].strip(),
                resnum=int(line[22:26]),
                icode=line[26],
                element=line[76:78].strip(),
                raw=line,
            ))

    return atoms

@dataclass
class Chain:
    key: str
    atom_lines: list = field(default_factory=list)
    residues: list = field(default_factory=list)
    seq: str = ""
    orig_chain_id: str = ""
    assigned_chain: str = ""
    split_reason: str = ""

def _add_line(chain, line):
    resnum = int(line[22:26])
    icode = line[26]
    resname = norm_resname(line[17:20].strip())
    chain.atom_lines.append(line)
    key = (resnum, icode)
    # residues가 비었거나, 직전 잔기와 (번호,삽입코드)가 다르면 → 새 잔기 시작
    if not chain.residues or (chain.residues[-1][0], chain.residues[-1][2]) != key:
        chain.residues.append((resnum, resname, icode))

def parse_pdb(path, gap = 1) :
    """PDB를 읽어 (헤더줄들, 사슬리스트, 경고리스트, 사용한분리방법) 반환."""
    header, warnings, raw = [], [], []
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "CRYST1" and not raw:
                header.append(line); continue
            if rec in ("ATOM", "HETATM", "TER"):
                raw.append(line)

    present = set()
    for line in raw:
        if line[:6].strip() in ("ATOM","HETATM") and line[17:20].strip() not in SKIP_RES:
            if line[21].strip():
                present.add(line[21])

    chains, method = [], ""
    def flush(c):
        if c.atom_lines: chains.append(c)

    if len(present) >= 2:
        method = "existing chain IDs"
        buckets = {}
        for line in raw:
            if line[:6].strip() not in ("ATOM","HETATM"): continue
            if line[17:20].strip() in SKIP_RES: continue
            cid = line[21]
            if cid not in buckets:
                buckets[cid] = Chain(key=cid, orig_chain_id=cid, split_reason="existing chain ID")
            _add_line(buckets[cid], line)
        chains = list(buckets.values())
    else:
        ter = sum(1 for line in raw if line[:6].strip() == "TER")
        if ter >= 1:
            method = "TER records"
            cur = Chain(key="0", split_reason="TER record")
            for line in raw:
                if line[:6].strip() == "TER":
                    flush(cur); cur = Chain(key=str(len(chains)), split_reason="TER record"); continue
                if line[:6].strip() not in ("ATOM","HETATM"): continue
                if line[17:20].strip() in SKIP_RES: continue
                _add_line(cur, line)
            flush(cur)
            if len(chains) < 2:
                chains = []
        if not chains:
            method = "residue-number discontinuity (FALLBACK)"
            warnings.append("No usable chain IDs or TER; chains inferred from residue-number jumps. VERIFY the split.")
            cur = Chain(key="0", split_reason="resnum discontinuity"); prev = None
            for line in raw:
                if line[:6].strip() not in ("ATOM","HETATM"): continue
                if line[17:20].strip() in SKIP_RES: continue
                num = int(line[22:26])
                if prev is not None and (num < prev or num - prev > gap + 1):
                    flush(cur); cur = Chain(key=str(len(chains)), split_reason="resnum discontinuity")
                _add_line(cur, line); prev = num
            flush(cur)

    for c in chains:
        c.seq = "".join(THREE_TO_ONE.get(rn, "X") for _, rn, _ in c.residues)
    if not chains:
        raise ValueError("No protein chains found.")
    return header, chains, warnings, method