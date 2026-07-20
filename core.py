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
                name=line[12:16].strip(),      # 원자 이름 컬럼
                resname=line[17:20].strip(),   # 잔기 이름 컬럼
                chain_id=line[21].strip(),     # 체인 ID 한 글자. 공백이면 "" 이 됨
                resnum=int(line[22:26]),       # 잔기 번호를 정수로 변환
                icode=line[26],                # 삽입 코드 (한 글자, 공백 그대로 둠)
                element=line[76:78].strip(),   # 원소 기호
                raw=line,                      # 원본 줄 보관 (커밋 6에서 재사용)
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
    """원자 줄 하나를 사슬에 추가. 잔기가 바뀌는 순간에만 residues에 새 항목을 넣는다."""
    resnum = int(line[22:26])
    icode = line[26]
    resname = norm_resname(line[17:20].strip())  # 잔기명 정규화 (norm_resname은 커밋 5에서 정의)
    chain.atom_lines.append(line)                # 원자 줄은 무조건 보관
    key = (resnum, icode)                        # 이 원자가 속한 잔기의 식별키
    # residues가 비었거나, 직전 잔기와 (번호,삽입코드)가 다르면 → 새 잔기 시작
    if not chain.residues or (chain.residues[-1][0], chain.residues[-1][2]) != key:
        chain.residues.append((resnum, resname, icode))