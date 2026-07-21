from dataclasses import dataclass, field
import string

THREE_TO_ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
                "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
                "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

NONSTD_MAP = {"HIE":"HIS","HID":"HIS","HIP":"HIS","HSD":"HIS","HSE":"HIS","HSP":"HIS",
              "CYX":"CYS","CYM":"CYS","ASH":"ASP","GLH":"GLU","LYN":"LYS","ARN":"ARG","MSE":"MET"}

SKIP_RES = {"HOH","WAT","TIP3","TIP","TIP4","SOL","NA","CL","SOD","CLA","POT",
            "K","MG","ZN","CA","MN","FE","CU","SO4","PO4","GOL","EDO"}

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

def norm_resname(resname):
    """비표준 잔기명이면 표준으로 바꾸고, 아니면 그대로 반환."""
    r = resname.strip()
    return NONSTD_MAP.get(r, r) 

def is_hydrogen(line):
    """이 원자 줄이 수소인가? 원소 컬럼을 우선 보고, 비어있으면 원자 이름으로 판단."""
    element = line[76:78].strip()
    atom_name = line[12:16].strip()
    if element:
        return element == "H"
    return atom_name.startswith("H") or (len(atom_name) > 1 and atom_name[0] in "123" and atom_name[1] == "H")

def write_clean(path, header, chains, strip_h=True):
    """정리된 사슬들을 PDB 형식으로 다시 써서 파일로 저장."""
    out = list(header)
    for c in chains:
        for line in c.atom_lines:
            if strip_h and is_hydrogen(line):
                continue
            resname = norm_resname(line[17:20].strip()).rjust(3)
            resnum = int(line[22:26])
            out.append(line[:17] + resname + " " + c.assigned_chain + "%4d" % resnum + line[26:])
        out.append("TER\n")
    out.append("END\n")
    with open(path, "w") as fh:
        fh.writelines(out)
    return path

def assign_chain_ids(chains, preferred=None):
    """각 사슬에 단일 글자 체인 ID를 배정. 원본이 전부 있고 서로 다르면 그대로, 아니면 A,B,C..."""
    origs = [c.orig_chain_id for c in chains]
    if all(origs) and len(set(origs)) == len(origs):
        for c in chains:
            c.assigned_chain = c.orig_chain_id
        return
    used = set(); letters = list(string.ascii_uppercase)
    for c in chains:
        cid = (preferred or {}).get(c.key) or next(l for l in letters if l not in used)
        c.assigned_chain = cid
        used.add(cid)