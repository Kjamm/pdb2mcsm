# pdb2mcsm

<img width="1248" height="704" alt="image" src="https://github.com/user-attachments/assets/fe188618-231b-4018-9cbd-34f75a5383e1" />

[English](https://github.com/Kjamm/pdb2mcsm/blob/main/README.md) &nbsp; &nbsp;[한국어](https://github.com/Kjamm/pdb2mcsm/blob/main/READMEKO.md)

Clean and prepare protein complex PDB files for [mCSM-PPI2](https://biosig.lab.uq.edu.au/mcsm_ppi2/) ΔΔG analysis.

MD simulation outputs and other "raw" structure files often can't be fed directly into mCSM-PPI2: chains may be merged or unlabeled, residues renumbered from 1, waters and ions cluttering the file, and histidines written in MD-specific names like `HIE`. `pdb2mcsm` fixes all of this automatically so the structure is accepted by the server.

## Try it in Colab (no install)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kjamm/pdb2mcsm/blob/main/pdb2mcsm_colab.ipynb)

Upload a PDB, get a cleaned `_ready.pdb` back, and optionally generate an alanine-scanning mutation list — all through form fields, no coding required.

## What it does

**Automatic cleanup**
- **Chain splitting** by priority: existing chain IDs → `TER` records → residue-number jumps (flagged when guessed).
- Handles **any number of chains** (dimer, trimer, antibody–antigen, ...).
- Removes **water, ions, and common crystallization additives** (HOH, NA, CL, SO4, GOL, ...).
- Strips **hydrogens** (mCSM re-adds them).
- Normalizes **MD residue names**: HIE/HID/HIP → HIS, CYX → CYS, MSE → MET, and more.
- Handles **phosphorylated residues**: SEP/S1P → SER, TPO/T1P → THR, PTR/Y1P → TYR, with phosphate atoms removed.
- Preserves **insertion codes** (e.g. 100A).
- Assigns clean chain IDs (keeps originals when valid).

**Optional**
- **Renumber** a chain by a fixed offset (e.g. shift file numbering to paper numbering).
- Generate an **alanine-scanning mutation list** in mCSM-PPI2 batch format.

**Safety** — the tool warns instead of silently guessing: ambiguous chain splits and single-chain inputs raise clear warnings so you can verify before trusting the output.

## Quick start (Python)

```python
import core

# 1. Read and split chains
header, chains, warnings, method = core.parse_pdb("my_complex.pdb")
core.assign_chain_ids(chains)

# 2. Review (important!)
print("split method:", method)
for row in core.summarize(chains):
    print(row)          # (chain, #residues, first, last, split_reason)
for w in warnings:
    print("WARNING:", w)

# 3. (optional) shift a chain's numbering: file 36 -> paper 98 means offset 62
core.renumber_chain(chains, "B", 62)

# 4. Write the cleaned file -> upload this to mCSM-PPI2
core.write_clean("my_complex_ready.pdb", header, chains)

# 5. (optional) alanine-scanning list for the batch box
print("\n".join(core.build_ala_scan(chains, "B", 95, 110)))
```

## Then in mCSM-PPI2

1. Open the mCSM-PPI2 batch/multiple-mutation page.
2. Upload `my_complex_ready.pdb`.
3. Select the chain to mutate.
4. Paste the alanine-scanning list.
5. Run → per-residue ΔΔG. Large negative values mark residues important for binding.

## Functions

| Function | Purpose |
|---|---|
| `parse_pdb(path)` | Read a PDB; split into chains; return `(header, chains, warnings, method)` |
| `assign_chain_ids(chains)` | Give each chain a clean A/B/C ID |
| `write_clean(path, header, chains)` | Write the cleaned PDB (H stripped, phosphate atoms removed) |
| `renumber_chain(chains, chain_id, offset)` | Shift a chain's residue numbers by `offset` |
| `build_ala_scan(chains, chain_id, start, end)` | Build an mCSM-PPI2 alanine-scanning list |
| `summarize(chains)` | Per-chain summary rows for review |
| `read_atoms(path)` | Low-level: parse ATOM/HETATM records into `Atom` objects |

## Limitations

`pdb2mcsm` solves file **structure and cleanup** problems. It does **not** handle:

- **Automatic renumbering** — you must know the offset yourself (`renumber_chain` is a fixed shift).
- **Other modified residues** beyond phosphorylation (acetylation, methylation, glycosylation, ...) — add them to `NONSTD_MAP` as needed.
- **Ligands / cofactors** (drugs, heme, ATP, metal clusters) — only water/ions are removed.
- **Structural defects** — missing atoms, broken backbone, atom clashes, or alternate locations (altloc) are not repaired.
- **Nucleic acids** — protein–DNA/RNA complexes are out of scope (mCSM-PPI2 is protein–protein only).

Some limits are **mCSM-PPI2's own**, not this tool's: phosphorylation *effects* are lost once a residue is converted to its standard form; protein–ligand and protein–nucleic-acid binding aren't supported; very large complexes may exceed server limits.

## Notes

- Default keeps **original residue numbering** — the numbers you feed mCSM must match what's in the cleaned file. Use `summarize()` or inspect `chains[i].residues` to confirm.
- Converting a phosphorylated residue to its standard form removes the phosphorylation from the model. If the phosphorylation itself matters for your binding question, interpret results accordingly or use a method that supports modified residues.

## License

MIT
