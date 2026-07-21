from dataclasses import dataclass, field
import string

THREE_TO_ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
                "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
                "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

NONSTD_MAP = {"HIE":"HIS","HID":"HIS","HIP":"HIS","HSD":"HIS","HSE":"HIS","HSP":"HIS",
              "CYX":"CYS","CYM":"CYS","ASH":"ASP","GLH":"GLU","LYN":"LYS","ARN":"ARG","MSE":"MET",
              "SEP":"SER","S1P":"SER",
              "TPO":"THR","T1P":"THR",
              "PTR":"TYR","Y1P":"TYR","PTY":"TYR"}

PHOSPHO_ATOMS = {"P","O1P","O2P","O3P","OP1","OP2","OP3","H1P","H2P","H3P","HOP2","HOP3"}

PHOSPHO_RES = {"SEP","S1P","TPO","T1P","PTR","Y1P","PTY"}

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

    if len(present) >= 1:
        method = "existing chain IDs"
        chains = []
        cur = None
        cur_cid = None
        prev_num = None
        for line in raw:
            rec = line[:6].strip()
            if rec == "TER":
                if cur: flush(cur); cur = None; prev_num = None
                continue
            if rec not in ("ATOM","HETATM"): continue
            if line[17:20].strip() in SKIP_RES: continue
            cid = line[21]
            num = int(line[22:26])
            if cur is None or cid != cur_cid or (prev_num is not None and num < prev_num):
                if cur: flush(cur)
                cur = Chain(key=str(len(chains)), orig_chain_id=cid, split_reason="existing chain ID")
                cur_cid = cid
            _add_line(cur, line)
            prev_num = num
        if cur: flush(cur)
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

    if len(chains) < 2:
        warnings.append(f"Only ONE chain detected (method: {method}). mCSM-PPI2 needs >=2 chains.")
    return header, chains, warnings, method

def norm_resname(resname):
    r = resname.strip()
    return NONSTD_MAP.get(r, r)

def is_hydrogen(line):
    element = line[76:78].strip()
    atom_name = line[12:16].strip()
    if element:
        return element == "H"
    return atom_name.startswith("H") or (len(atom_name) > 1 and atom_name[0] in "123" and atom_name[1] == "H")

def write_clean(path, header, chains, strip_h=True):
    out = list(header)
    for c in chains:
        for line in c.atom_lines:
            if strip_h and is_hydrogen(line):
                continue
            if line[17:20].strip() in PHOSPHO_RES and line[12:16].strip() in PHOSPHO_ATOMS:
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

def build_ala_scan(chains, chain_id, start, end):
    target = None
    for c in chains:
        if c.assigned_chain == chain_id:
            target = c; break
    if target is None:
        raise ValueError(f"chain {chain_id} not found")
    muts = []
    for resnum, resname, icode in target.residues:
        if not (start <= resnum <= end):
            continue
        one = THREE_TO_ONE.get(resname, "X")
        if one == "A":
            muts.append(f"{chain_id} A{resnum}G")
        elif one == "G":
            muts.append(f"{chain_id} G{resnum}A")
        elif one != "X":
            muts.append(f"{chain_id} {one}{resnum}A")
    return muts

def summarize(chains):
    rows = []
    for c in chains:
        first = f"{c.residues[0][1]}{c.residues[0][0]}"
        last = f"{c.residues[-1][1]}{c.residues[-1][0]}"
        rows.append((c.assigned_chain or c.key, len(c.residues), first, last, c.split_reason))
    return rows

def renumber_chain(chains, chain_id, offset):
    """지정한 사슬의 모든 잔기 번호에 offset을 더한다(빼려면 음수). 
    파일번호 → 논문번호처럼 번호 체계를 옮길 때 사용.
    예: 파일에서 36번인 잔기를 98번으로 만들려면 offset=62.

    chains 안의 원자 줄(atom_lines)과 잔기 목록(residues)을 모두 갱신하므로,
    이 함수 호출 뒤 write_clean을 하면 새 번호로 저장된다."""
    target = None
    for c in chains:
        if c.assigned_chain == chain_id or c.orig_chain_id == chain_id:
            target = c; break
    if target is None:
        raise ValueError(f"chain {chain_id} not found")

    new_lines = []
    for line in target.atom_lines:
        old_num = int(line[22:26])
        new_num = old_num + offset
        new_lines.append(line[:22] + "%4d" % new_num + line[26:])
    target.atom_lines = new_lines

    target.residues = [(num + offset, resname, icode)
                       for (num, resname, icode) in target.residues]
    return target