# Phasing work: notes and run provenance

Kept under version control because they existed only in the gitignored
`test_data/` and `test_data2/` directories.

| file | what it is |
|---|---|
| `phasing-roadmap.md` | the running record of the phasing analysis: what was tried, what was measured, what is open |
| `joint-phasing-brief.md` | design brief for joint cross-cell phasing (C2) |
| `runs/S1/command*.sh` | the exact command for each S1 run; `commandN.sh` wrote `hscn_N.csv.gz`, except `command.sh`, which wrote `hscn_relabel.csv.gz` |
| `runs/S2/command*.sh` | the same for S2 (`command.sh` -> r1, `command_c2.sh` -> c2, `command3.sh` -> r3, `command4.sh` -> r4) |

`bench/history.tsv` maps these runs to named configurations and records which
fixes each had - the same flags meant different things at different script
versions.

These are **copies**. The command scripts must stay in `test_data/` and
`test_data2/`, since they run from there with relative paths. The roadmap and
brief were last edited in `test_data/`; edit the copies here from now on so the
two do not drift apart.

They name the samples, so they should not be pushed to this public repository
without deciding that is acceptable.
