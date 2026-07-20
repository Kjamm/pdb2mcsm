from dataclasses import dataclass

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

a = read_atoms("/Users/jaeminkim/pdb2mcsm/pdb2mcsm/tests/fixtures/mini_raw.pdb")
print(len(a))                    # 원자 총 개수
print(a[0].resname, a[0].resnum) # 첫 잔기 이름/번호가 파일과 일치하는지