# DANDI 000003 Targeted Download Priority

Date built: 2026-06-24

## Rationale

The default multi-session download plan ranks missing NWB files by size. This targeted plan instead asks which single additional file is most likely to improve the dentate-granule physiology validation layer per GB.

Scoring favors: same-subject sessions from animals that already produced labeled granule units, adjacent recording days, small file size, and a smaller bonus for new-subject breadth. The score is heuristic and should guide downloads, not replace source-level unit validation.

## Recommended Next File

- Session: `YutaMouse55-160911`
- Subject: `YutaMouse55`
- Asset: `5124caca-1174-44e2-8529-f306bd383a2e`
- Size: 9.01 GB
- Track: `yield_first`
- Reason: same subject has 16 local granule units; 2-day gap from local session; 9.01 GB

This is preferred over the smallest missing file because it is an adjacent follow-up from a subject that already yielded labeled granule units locally.

## Top Candidates

| Rank | Session | Subject | Track | Size GB | Score/GB | Reason |
|---:|---|---|---|---:|---:|---|
| 1 | `YutaMouse55-160911` | `YutaMouse55` | `yield_first` | 9.01 | 1.695 | same subject has 16 local granule units; 2-day gap from local session; 9.01 GB |
| 2 | `YutaMouse42-151103` | `YutaMouse42` | `yield_first` | 7.77 | 1.684 | same subject has 2 local granule units; 1-day gap from local session; 7.77 GB |
| 3 | `YutaMouse37-150610` | `YutaMouse37` | `yield_first` | 10.09 | 1.559 | same subject has 5 local granule units; 1-day gap from local session; 10.09 GB |
| 4 | `YutaMouse55b160907` | `YutaMouse55` | `yield_first` | 11.78 | 1.323 | same subject has 16 local granule units; 0-day gap from local session; 11.78 GB |
| 5 | `YutaMouse55-160910` | `YutaMouse55` | `yield_first` | 11.78 | 1.323 | same subject has 16 local granule units; 1-day gap from local session; 11.78 GB |
| 6 | `YutaMouse55-160906` | `YutaMouse55` | `yield_first` | 12.10 | 1.282 | same subject has 16 local granule units; 1-day gap from local session; 12.10 GB |
| 7 | `YutaMouse42-151117` | `YutaMouse42` | `yield_first` | 8.38 | 1.263 | same subject has 2 local granule units; 15-day gap from local session; 8.38 GB |
| 8 | `YutaMouse42-151114` | `YutaMouse42` | `yield_first` | 9.26 | 1.195 | same subject has 2 local granule units; 12-day gap from local session; 9.26 GB |
| 9 | `YutaMouse41-150901` | `YutaMouse41` | `yield_first` | 11.76 | 1.129 | same subject has 3 local granule units; 3-day gap from local session; 11.76 GB |
| 10 | `YutaMouse55-160902` | `YutaMouse55` | `yield_first` | 12.54 | 1.110 | same subject has 16 local granule units; 5-day gap from local session; 12.54 GB |
| 11 | `YutaMouse37-150612` | `YutaMouse37` | `yield_first` | 13.89 | 1.024 | same subject has 5 local granule units; 3-day gap from local session; 13.89 GB |
| 12 | `YutaMouse41-150830` | `YutaMouse41` | `yield_first` | 13.51 | 1.007 | same subject has 3 local granule units; 1-day gap from local session; 13.51 GB |

## Output

- Full priority table: `Project/results/dandi_000003_targeted_download_priority.tsv`
