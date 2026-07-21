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

    if len(present) >= 1:
        # 경로 A: 파일에 chain ID가 있으면 존중.
        # 단, "같은 ID인데 TER/번호리셋으로 실제로는 여러 사슬"인 경우까지 잡기 위해
        # ID별 1차 그룹핑 후, 각 그룹 안에서 TER 또는 번호 되돌아감으로 2차 분리한다.
        method = "existing chain IDs"
        chains = []
        cur = None            # 현재 쌓는 중인 사슬
        cur_cid = None        # 현재 원본 chain id
        prev_num = None       # 직전 잔기 번호 (같은 그룹 내 리셋 감지용)
        for line in raw:
            rec = line[:6].strip()
            if rec == "TER":              # TER는 무조건 사슬 경계
                if cur: flush(cur); cur = None; prev_num = None
                continue
            if rec not in ("ATOM","HETATM"): continue
            if line[17:20].strip() in SKIP_RES: continue
            cid = line[21]
            num = int(line[22:26])
            # 새 사슬을 시작해야 하는가?
            #  - 아직 사슬이 없음, 또는
            #  - chain id가 바뀜, 또는
            #  - 같은 id인데 번호가 되돌아감(리셋) → 실제로는 다른 분자
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

def build_ala_scan(chains, chain_id, start, end):
    """지정한 사슬의 [start, end] 구간 잔기를 alanine scanning 형식 목록으로 만든다."""
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
    """사람이 확인하기 좋게 사슬별 요약 행을 만든다: (체인, 잔기수, 첫잔기, 끝잔기, 분리근거)."""
    rows = []
    for c in chains:
        first = f"{c.residues[0][1]}{c.residues[0][0]}"
        last = f"{c.residues[-1][1]}{c.residues[-1][0]}"
        rows.append((c.assigned_chain or c.key, len(c.residues), first, last, c.split_reason))
    return rows