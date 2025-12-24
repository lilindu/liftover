# pombe_liftover

A static, client‑side liftover tool that converts genomic coordinates between the PomBase reference genome and the Leupold consensus genome. The app runs entirely in the browser and is published via GitHub Pages.

- Live site: `https://lilindu.github.io/liftover/`
- Local Pages copy: `http://localhost:8001/docs/`
- Local dev app: `http://localhost:8000/web/`

## Repository Layout

- `docs/index.html` — Web UI and styles for the GitHub Pages site
- `docs/app.js` — Liftover logic, UI wiring, status updates, version badge
- `docs/A_to_B.blocks.tsv` — Blocks mapping A→B (PomBase → Leupold)
- `docs/B_to_A.blocks.tsv` — Blocks mapping B→A (Leupold → PomBase)
- `liftover_tool/web/index.html` — Local dev HTML
- `liftover_tool/web/app.js` — Local dev JS

## Blocks File Generation

Goal: Produce high‑confidence, contiguous A→B and B→A blocks that map the PomBase reference to the Leupold consensus with strand and contig annotations suitable for liftover.

Inputs:
- A (source): PomBase reference chromosomes I/II/III (FASTA)
- B (target): Leupold consensus chromosomes I/II/III (FASTA)

**Run Alignment**
- Command:
  - `minimap2 -cx asm5 -c --secondary=no Leupold_consensus_genome_250919/Lcon_v250825_genome.fa Schizosaccharomyces_pombe_chromosome_I.fa Schizosaccharomyces_pombe_chromosome_II.fa Schizosaccharomyces_pombe_chromosome_III.fa > A_to_B.paf`

**Generate Mapping Blocks**
- Command:
  - `python3 cli/paf_to_blocks.py A_to_B.paf -o A_to_B.blocks.tsv`
- Block schema:
  - `contigA  startA  endA  contigB  startB  endB  strand  mapq` (0‑based, half‑open)

Pipeline notes:
- Use minimap2 `-x asm5` for high‑identity assemblies; `--secondary=no` filters secondary matches.
- Merge adjacent A→B hits on the same `contigB` and `strand` while enforcing A‑side non‑overlap.
- Confirm homogeneity (target contig and strand) within each block.
- Sort blocks per `contigA` by `startA`.

## Coordinate Systems

### Coordinate Systems at a Glance

| Item | Blocks TSV (0-based, half-open) | UI input/output (1-based, inclusive) | Conversion |
|---|---|---|---|
| Point | `pos0` | `pos1` | `pos1 = pos0 + 1` |
| Interval start | `start0` | `start1` | `start1 = start0 + 1` |
| Interval end | `end0` (exclusive) | `end1` (inclusive) | `end1 = end0` |
| Interval length | `end0 - start0` | `end1 - start1 + 1` | equal |

### Examples

- Block row (TSV):
  - `contigA I  startA0 100  endA0 200  contigB …  startB0 5000  endB0 5100  strand +`
  - Visible span in UI terms: `I:101–200` maps to `…:5001–5100`.
- UI input `I:101–120` corresponds to TSV interval `[100,120)` on A; if mapped B endpoints in TSV are `5001–5020`, the UI shows `…:5001–5020`.
- Strand handling: TSV stores endpoints as 0‑based, half‑open; the strand affects how `mapPoint` computes B positions, but UI rendering remains 1‑based, inclusive.

## Liftover Algorithm

- Parse and index blocks
  - `parseTSV(text)` builds `{ contigA: Block[] }` sorted by `startA`.
- Block search
  - `findBlock(blocks, pos)` returns the index of the block containing `pos`.
- Point mapping
  - `mapPoint(block, posA)` maps a single A position to B respecting block strand.
- Single‑block interval
  - `mapInterval(block, startA, endA)` → `[contigB, startB, endB, strand]`.
- Stitching across blocks
  - `stitchInterval(blocks, startA, endA)`:
    - Fails as `UNMAPPED` if endpoints lie outside blocks.
    - Enforces target `contigB` and `strand` consistency.
    - Maps stitched ends via `mapPoint` and normalizes B `[startB,endB]`.
- Gap and overlap reporting
  - `collectGaps(blocks, startIdx, endIdx)` emits:
    - `GAP_A:N@A:start→end`, `GAP_B:N@B:start→end`, `OVERLAP_B:N@B:start→end`
    - `CONTIG_CHANGE`, `STRAND_CHANGE` when encountered
  - `formatGaps(gaps)` renders compact tokens.
- Execution
  - AB liftover iterates inputs; status = `OK` or `SPANNING_BLOCKS`, `UNMAPPED` on failure.
  - BA liftover mirrors AB.
- Round‑Trip validation
  - A→B then B→A; status `PASS` only if mapped‑back interval equals original, else `FAIL`.

## Status Semantics

- `OK` — Interval maps within a single block
- `SPANNING_BLOCKS` — Interval spans multiple blocks; see “gaps”
- `UNMAPPED` — One or both endpoints do not fall in any block
- Round‑Trip: `PASS` if B→A equals original; otherwise `FAIL`

## Web UI

- Title/subtitle: pombe_liftover, with a version badge (e.g., `v0.6 (2025-12-25)`).
- Status chips: AB/BA loaded blocks.
- Input area: textarea accepting `CHR:POS` or `CHR:START-END` (1‑based, inclusive).
- Buttons: AB, BA, Round‑Trip, Download (unified primary style).
- Explanation box: A=source; B=target; status and gaps syntax.
- Dynamic direction:
  - AB: `Current direction: A = PomBase → B = Leupold`
  - BA: `Current direction: A = Leupold → B = PomBase`
  - Round‑Trip: `Current direction: A = PomBase → B = Leupold → A = PomBase`
- Results table: sticky header, zebra rows, right‑aligned numerics, wrapped “gaps”.

## Versioning & Cache Busting

- `APP_VERSION` is a literal string (e.g., `v0.6 (2025-12-25)`) set in `app.js`.
- `index.html` loads `app.js` with `?v=autoload_fix_<n>` to force fresh loads on deployment.

## Local Development

- Serve `docs/`: `python3 -m http.server 8001` then open `http://localhost:8001/docs/`.
- Serve local app: `python3 -m http.server 8000` in `liftover_tool/` then open `http://localhost:8000/web/`.
- Testing:
  - AB/BA liftover produces `OK` / `SPANNING_BLOCKS` with meaningful `gaps`.
  - Round‑Trip produces `PASS/FAIL` with `gapsAB/gapsBA`.

## Deployment (GitHub Pages)

- Settings: Repository → Settings → Pages → Source: `main`, Folder: `/docs`.
- Push changes to `main`; Pages rebuild in ~1–2 minutes.
- Hard refresh `https://lilindu.github.io/liftover/`.

## Design Rationale

- Client‑only tool: easy hosting and reproducibility.
- Non‑overlapping A blocks: simplifies search and avoids ambiguity.
- Reporting gaps/overlaps on B: makes local alignment artifacts explicit.

## Extensibility

- Stitch tolerance slider to treat tiny gaps as contiguous.
- Table filters/sorting and copy helpers.
- URL state for sharing input sessions.
- Accessibility: keyboard shortcuts; high‑contrast theme toggle.

## Known Limitations

- Overlapping A blocks are not supported by design.
- If contig/strand transitions are encountered mid‑stitch, liftover fails to preserve correctness.
- Inputs must be 1‑based, inclusive; malformed lines are rejected.

## Code References

- Parse TSV: `liftover/docs/app.js:18–40`
- Binary search: `liftover/docs/app.js:42–53`
- Point mapping: `liftover/docs/app.js:55–66`
- Stitch interval: `liftover/docs/app.js:151–180`
- Collect gaps: `liftover/docs/app.js:182–220`
- Format gaps: `liftover/docs/app.js:222–233`
- AB liftover: `liftover/docs/app.js:264–295`
- BA liftover: `liftover/docs/app.js:297–328`
- Round‑Trip: `liftover/docs/app.js:330–369`
- Direction note update: `liftover/docs/app.js:427–436`
- Render TSV: `liftover/docs/app.js:235–262`
