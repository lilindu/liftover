// Global state
let blocksAB = null;
let blocksBA = null;
const APP_VERSION = 'v0.4 (2025-11-22)';

// Log to console and to the autoload debug panel
function debugLog(message) {
  const div = document.getElementById('autoloadDebug');
  if (div) {
    const ts = new Date().toLocaleTimeString();
    div.innerHTML += `<br>[${ts}] ${message}`;
    div.scrollTop = div.scrollHeight;
  }
  console.log('[AUTOLOAD]', message);
}

// Parse blocks TSV into per-contig sorted arrays
function parseTSV(text) {
  const lines = text.trim().split(/\r?\n/);
  const by = {};
  for (let i = 1; i < lines.length; i++) {
    const parts = lines[i].split("\t");
    if (parts.length < 8) continue;
    const b = {
      contigA: parts[0],
      startA: parseInt(parts[1]),
      endA: parseInt(parts[2]),
      contigB: parts[3],
      startB: parseInt(parts[4]),
      endB: parseInt(parts[5]),
      strand: parts[6],
      mapq: parseInt(parts[7])
    };
    if (!by[b.contigA]) by[b.contigA] = [];
    by[b.contigA].push(b);
  }
  for (const c in by) by[c].sort((a,b)=>a.startA-b.startA);
  return by;
}

// Binary search to find containing block for a position
function findBlock(blocks, pos) {
  let lo = 0, hi = blocks.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    const b = blocks[mid];
    if (pos < b.startA) hi = mid;
    else if (pos >= b.endA) lo = mid + 1;
    else return mid;
  }
  return null;
}

// Map position within a block from A to B according to strand
function mapPoint(block, posA) {
  const offset = posA - block.startA;
  if (block.strand === '+') {
    const posB = block.startB + offset;
    return [block.contigB, posB, block.strand];
  } else {
    const length = block.endA - block.startA;
    const posB = block.startB + (length - 1 - offset);
    return [block.contigB, posB, block.strand];
  }
}

// Fetch default PomBase→Leupold blocks from server with cache-busting
async function loadDefaultAB() {
  const statusAB = document.getElementById('statusAB');
  try {
    debugLog('Fetching PomBase→Leupold blocks...');
    const res = await fetch('A_to_B.blocks.tsv?v=' + Date.now());
    if (!res.ok) throw new Error('HTTP ' + res.status);
    const txt = await res.text();
    blocksAB = parseTSV(txt);
    const count = Object.keys(blocksAB).reduce((s, c) => s + blocksAB[c].length, 0);
    if (statusAB) statusAB.textContent = 'PomBase→Leupold: blocks file loaded (' + count + ' blocks)';
    debugLog('PomBase→Leupold blocks file loaded: ' + count + ' blocks');
  } catch (e) {
    if (statusAB) statusAB.textContent = 'PomBase→Leupold: failed (' + e.message + ')';
    debugLog('PomBase→Leupold failed: ' + e.message);
  }
}

// Fetch default Leupold→PomBase blocks from server with cache-busting
async function loadDefaultBA() {
  const statusBA = document.getElementById('statusBA');
  try {
    debugLog('Fetching Leupold→PomBase blocks...');
    const res = await fetch('B_to_A.blocks.tsv?v=' + Date.now());
    if (!res.ok) throw new Error('HTTP ' + res.status);
    const txt = await res.text();
    blocksBA = parseTSV(txt);
    const count = Object.keys(blocksBA).reduce((s, c) => s + blocksBA[c].length, 0);
    if (statusBA) statusBA.textContent = 'Leupold→PomBase: blocks file loaded (' + count + ' blocks)';
    debugLog('Leupold→PomBase blocks file loaded: ' + count + ' blocks');
  } catch (e) {
    if (statusBA) statusBA.textContent = 'Leupold→PomBase: failed (' + e.message + ')';
    debugLog('Leupold→PomBase failed: ' + e.message);
  }
}


// Parse one input line which can be a point (CHR:POS) or interval (CHR:START-END)
// Returns { contig, start, end, isInterval } or null on bad input
function parseInputLine(line) {
  const parts = line.split(':');
  if (parts.length !== 2) return null;
  const contig = parts[0];
  const coords = parts[1];
  if (/^\d+$/.test(coords)) {
    const pos = parseInt(coords);
    return { contig, start: pos, end: pos, isInterval: false };
  }
  const m = coords.match(/^(\d+)-(\d+)$/);
  if (!m) return null;
  const a = parseInt(m[1]);
  const b = parseInt(m[2]);
  const start = Math.min(a, b);
  const end = Math.max(a, b);
  return { contig, start, end, isInterval: true };
}

// Map an interval entirely contained within a single block; returns [contigB, startB, endB, strand]
function mapInterval(block, startA, endA) {
  const [, sB] = mapPoint(block, startA);
  const [, eB] = mapPoint(block, endA);
  const startB = Math.min(sB, eB);
  const endB = Math.max(sB, eB);
  return [block.contigB, startB, endB, block.strand];
}

// Attempt to stitch an interval across multiple blocks under strict constraints.
// Constraints: same target contig and strand across blocks; no gaps in A; monotonic and gap-free in B within tolerance.
// Returns { ok, reason, contigB, startB, endB, strand }
const STITCH_TOLERANCE = 0;
function stitchInterval(blocks, startA, endA, tolerance = STITCH_TOLERANCE) {
  if (startA > endA) { const t = startA; startA = endA; endA = t; }
  const idxStart = findBlock(blocks, startA);
  const idxEnd = findBlock(blocks, endA);
  if (idxStart === null || idxEnd === null) return { ok: false, reason: 'UNMAPPED' };

  const first = blocks[idxStart];
  const last = blocks[idxEnd];
  if (idxStart === idxEnd) {
    const [contigB, sB, eB, strand] = mapInterval(first, startA, endA);
    return { ok: true, contigB, startB: sB, endB: eB, strand };
  }

  const targetContig = first.contigB;
  const targetStrand = first.strand;

  for (let i = idxStart; i < idxEnd; i++) {
    const cur = blocks[i];
    const next = blocks[i + 1];
    if (cur.contigB !== targetContig) return { ok: false, reason: 'CONTIG' };
    if (cur.strand !== targetStrand) return { ok: false, reason: 'STRAND' };
    // Allow gaps; they will be reported separately via collectGaps
  }

  const sB = mapPoint(first, startA)[1];
  const eB = mapPoint(last, endA)[1];
  const startB = Math.min(sB, eB);
  const endB = Math.max(sB, eB);
  return { ok: true, contigB: targetContig, startB, endB, strand: targetStrand };
}

function collectGaps(blocks, startIdx, endIdx) {
  const gaps = [];
  const first = blocks[startIdx];
  const targetStrand = first.strand;
  for (let i = startIdx; i < endIdx; i++) {
    const cur = blocks[i];
    const next = blocks[i + 1];
    if (cur.contigB !== next.contigB) {
      gaps.push({ type: 'CONTIG_CHANGE', contigPrev: cur.contigB, contigNext: next.contigB });
    }
    if (cur.strand !== next.strand) {
      gaps.push({ type: 'STRAND_CHANGE', strandPrev: cur.strand, strandNext: next.strand });
    }
    const gapA = Math.max(0, next.startA - cur.endA);
    if (gapA > 0) {
      gaps.push({ type: 'GAP_A', sizeA: gapA, aStart: cur.endA + 1, aEnd: next.startA - 1 });
    }
    const bEnd = mapPoint(cur, cur.endA - 1)[1];
    const bStartNext = mapPoint(next, next.startA)[1];
    if (targetStrand === '+') {
      const gapB = Math.max(0, bStartNext - bEnd - 1);
      const ovlB = Math.max(0, bEnd - bStartNext + 1);
      if (gapB > 0) gaps.push({ type: 'GAP_B', sizeB: gapB, bStart: bEnd + 1, bEnd: bStartNext - 1 });
      if (ovlB > 0) gaps.push({ type: 'OVERLAP_B', sizeB: ovlB, bStart: bStartNext, bEnd: bEnd });
    } else {
      const gapB = Math.max(0, bEnd - bStartNext - 1);
      const ovlB = Math.max(0, bStartNext - bEnd + 1);
      if (gapB > 0) {
        const gs = Math.min(bEnd - 1, bStartNext + 1);
        const ge = Math.max(bEnd - 1, bStartNext + 1);
        gaps.push({ type: 'GAP_B', sizeB: gapB, bStart: gs, bEnd: ge });
      }
      if (ovlB > 0) {
        const os = Math.min(bStartNext, bEnd);
        const oe = Math.max(bStartNext, bEnd);
        gaps.push({ type: 'OVERLAP_B', sizeB: ovlB, bStart: os, bEnd: oe });
      }
    }
  }
  return gaps;
}

function formatGaps(gaps) {
  if (!gaps || gaps.length === 0) return '';
  return gaps.map(g => {
    if (g.type === 'GAP_A') return `GAP_A:${g.sizeA}@A:${g.aStart}→${g.aEnd}`;
    if (g.type === 'GAP_B') return `GAP_B:${g.sizeB}@B:${g.bStart}→${g.bEnd}`;
    if (g.type === 'OVERLAP_B') return `OVERLAP_B:${g.sizeB}@B:${g.bStart}→${g.bEnd}`;
    if (g.type === 'CONTIG_CHANGE') return `CONTIG_CHANGE@B:${g.contigPrev}→${g.contigNext}`;
    if (g.type === 'STRAND_CHANGE') return `STRAND_CHANGE@B:${g.strandPrev}→${g.strandNext}`;
    return g.type;
  }).join('; ');
}

// Render a TSV string into the results table while keeping the raw text for download
function renderTSVToTable(tsv) {
  const head = document.getElementById('outHead');
  const body = document.getElementById('outBody');
  if (!head || !body) return;
  head.innerHTML = '';
  body.innerHTML = '';
  const lines = tsv.trim().split(/\r?\n/);
  if (lines.length === 0) return;
  const headers = lines[0].split('\t');
  const trh = document.createElement('tr');
  for (const h of headers) {
    const th = document.createElement('th');
    th.textContent = h;
    trh.appendChild(th);
  }
  head.appendChild(trh);
  for (let i = 1; i < lines.length; i++) {
    const cols = lines[i].split('\t');
    const tr = document.createElement('tr');
    for (let j = 0; j < headers.length; j++) {
      const td = document.createElement('td');
      td.textContent = cols[j] || '';
      tr.appendChild(td);
    }
    body.appendChild(tr);
  }
}

// Run liftover PomBase→Leupold for input lines (supports points and intervals)
function runLiftoverAB() {
  const out = [];
  out.push('contigA\tstartA\tendA\tcontigB\tstartB\tendB\tstrand\tstatus\tgaps');
  if (!blocksAB) {
    out.push('\t\t\t\t\t\t\tNO_BLOCKS\t');
    document.getElementById('out').textContent = out.join('\n');
    return;
  }
  const lines = document.getElementById('coords').value.trim().split(/\r?\n/);
  for (const s of lines) {
    if (!s) continue;
    const parsed = parseInputLine(s);
    if (!parsed) { out.push('\t\t\t\t\t\t\tBAD_INPUT\t'); continue; }
    const { contig, start, end } = parsed;
    const blocks = blocksAB[contig];
    if (!blocks) { out.push(`${contig}\t${start}\t${end}\t\t\t\t\tNO_CONTIG`); continue; }
    const attempt = stitchInterval(blocks, start, end);
    if (!attempt.ok) {
      const statusFail = (attempt.reason === 'UNMAPPED') ? 'UNMAPPED' : `STITCH_FAILED_${attempt.reason}`;
      out.push(`${contig}\t${start}\t${end}\t\t\t\t\t${statusFail}\t`);
      continue;
    }
    const sameBlock = (findBlock(blocks, start) === findBlock(blocks, end));
    const gaps = sameBlock ? [] : collectGaps(blocks, findBlock(blocks, start), findBlock(blocks, end));
    const status = sameBlock ? 'OK' : 'SPANNING_BLOCKS';
    out.push(`${contig}\t${start}\t${end}\t${attempt.contigB}\t${attempt.startB}\t${attempt.endB}\t${attempt.strand}\t${status}\t${formatGaps(gaps)}`);
  }
  document.getElementById('out').textContent = out.join('\n');
  renderTSVToTable(document.getElementById('out').textContent);
}

// Run liftover Leupold→PomBase for input lines (supports points and intervals)
function runLiftoverBA() {
  const out = [];
  out.push('contigA\tstartA\tendA\tcontigB\tstartB\tendB\tstrand\tstatus\tgaps');
  if (!blocksBA) {
    out.push('\t\t\t\t\t\t\tNO_BLOCKS\t');
    document.getElementById('out').textContent = out.join('\n');
    return;
  }
  const lines = document.getElementById('coords').value.trim().split(/\r?\n/);
  for (const s of lines) {
    if (!s) continue;
    const parsed = parseInputLine(s);
    if (!parsed) { out.push('\t\t\t\t\t\t\tBAD_INPUT\t'); continue; }
    const { contig, start, end } = parsed;
    const blocks = blocksBA[contig];
    if (!blocks) { out.push(`${contig}\t${start}\t${end}\t\t\t\t\tNO_CONTIG`); continue; }
    const attempt = stitchInterval(blocks, start, end);
    if (!attempt.ok) {
      const statusFail = (attempt.reason === 'UNMAPPED') ? 'UNMAPPED' : `STITCH_FAILED_${attempt.reason}`;
      out.push(`${contig}\t${start}\t${end}\t\t\t\t\t${statusFail}\t`);
      continue;
    }
    const sameBlock = (findBlock(blocks, start) === findBlock(blocks, end));
    const gaps = sameBlock ? [] : collectGaps(blocks, findBlock(blocks, start), findBlock(blocks, end));
    const status = sameBlock ? 'OK' : 'SPANNING_BLOCKS';
    out.push(`${contig}\t${start}\t${end}\t${attempt.contigB}\t${attempt.startB}\t${attempt.endB}\t${attempt.strand}\t${status}\t${formatGaps(gaps)}`);
  }
  document.getElementById('out').textContent = out.join('\n');
  renderTSVToTable(document.getElementById('out').textContent);
}

// Run round-trip PomBase→Leupold→PomBase for input lines (supports points and intervals)
function runRoundtrip() {
  const out = [];
  out.push('contigA\tstartA\tendA\tcontigB\tstartB\tendB\tcontigA2\tstartA2\tendA2\tstatus\tgapsAB\tgapsBA');
  if (!blocksAB || !blocksBA) {
    out.push('\t\t\t\t\t\t\t\t\tNO_BLOCKS\t\t');
    document.getElementById('out').textContent = out.join('\n');
    return;
  }
  const lines = document.getElementById('coords').value.trim().split(/\r?\n/);
  for (const s of lines) {
    if (!s) continue;
    const parsed = parseInputLine(s);
    if (!parsed) { out.push('\t\t\t\t\t\t\t\t\tBAD_INPUT\t\t'); continue; }
    const { contig, start, end } = parsed;

    // Stage 1: PomBase→Leupold stitch
    const blocks1 = blocksAB[contig];
    if (!blocks1) { out.push(`${contig}\t${start}\t${end}\t\t\t\t\t\t\tNO_CONTIG`); continue; }
    const attempt1 = stitchInterval(blocks1, start, end);
    if (!attempt1.ok) { out.push(`${contig}\t${start}\t${end}\t\t\t\t\t\t\tSTITCH_FAILED_${attempt1.reason}\t\t`); continue; }
    const sameAB = (findBlock(blocks1, start) === findBlock(blocks1, end));
    const gapsAB = sameAB ? [] : collectGaps(blocks1, findBlock(blocks1, start), findBlock(blocks1, end));

    const blocks2 = blocksBA[attempt1.contigB];
    if (!blocks2) { out.push(`${contig}\t${start}\t${end}\t${attempt1.contigB}\t${attempt1.startB}\t${attempt1.endB}\t\t\t\tNO_CONTIG_BA`); continue; }
    const attempt2 = stitchInterval(blocks2, attempt1.startB, attempt1.endB);
    if (!attempt2.ok) { const status2 = (attempt2.reason === 'UNMAPPED') ? 'UNMAPPED_BA' : `STITCH_FAILED_BA_${attempt2.reason}`; out.push(`${contig}\t${start}\t${end}\t${attempt1.contigB}\t${attempt1.startB}\t${attempt1.endB}\t\t\t\t${status2}\t${formatGaps(gapsAB)}\t`); continue; }
    const sameBA = (findBlock(blocks2, attempt1.startB) === findBlock(blocks2, attempt1.endB));
    const gapsBA = sameBA ? [] : collectGaps(blocks2, findBlock(blocks2, attempt1.startB), findBlock(blocks2, attempt1.endB));

    const startA2 = Math.min(attempt2.startB, attempt2.endB);
    const endA2 = Math.max(attempt2.startB, attempt2.endB);
    const status = (attempt2.contigB === contig && startA2 === start && endA2 === end) ? 'PASS' : 'FAIL';
    out.push(`${contig}\t${start}\t${end}\t${attempt1.contigB}\t${attempt1.startB}\t${attempt1.endB}\t${attempt2.contigB}\t${startA2}\t${endA2}\t${status}\t${formatGaps(gapsAB)}\t${formatGaps(gapsBA)}`);
  }
  document.getElementById('out').textContent = out.join('\n');
  renderTSVToTable(document.getElementById('out').textContent);
}

// Download results TSV
function downloadResults() {
  const text = document.getElementById('out').textContent;
  const blob = new Blob([text], { type: 'text/tab-separated-values' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = 'liftover_results.tsv'; a.click();
  URL.revokeObjectURL(url);
}

// Initialize: attach handlers and autoload blocks
async function init() {
  debugLog('=== Starting autoload initialization ===');
  debugLog('DOM readyState: ' + document.readyState);

  function writeVersion() {
    const vb = document.getElementById('versionBadge');
    if (vb) vb.textContent = APP_VERSION;
  }
  writeVersion();
  const vt = document.getElementById('versionText');
  if (vt) vt.textContent = APP_VERSION;

  const sAB = document.getElementById('statusAB');
  const sBA = document.getElementById('statusBA');
  if (sAB) sAB.textContent = 'PomBase→Leupold: autoload starting';
  if (sBA) sBA.textContent = 'Leupold→PomBase: autoload starting';

  // Attach UI handlers
  const el = (id) => document.getElementById(id);
  el('loadDefaultAB') && el('loadDefaultAB').addEventListener('click', loadDefaultAB);
  el('loadDefaultBA') && el('loadDefaultBA').addEventListener('click', loadDefaultBA);
  el('runAB') && el('runAB').addEventListener('click', runLiftoverAB);
  el('runBA') && el('runBA').addEventListener('click', runLiftoverBA);
  el('runRoundtrip') && el('runRoundtrip').addEventListener('click', runRoundtrip);
  el('download') && el('download').addEventListener('click', downloadResults);
  el('forceAutoload') && el('forceAutoload').addEventListener('click', init);

  try {
    debugLog('Calling loadDefaultAB...');
    await loadDefaultAB();
    debugLog('loadDefaultAB completed');
  } catch (e) {
    debugLog('loadDefaultAB failed: ' + e.message);
  }
  try {
    debugLog('Calling loadDefaultBA...');
    await loadDefaultBA();
    debugLog('loadDefaultBA completed');
  } catch (e) {
    debugLog('loadDefaultBA failed: ' + e.message);
  }

  debugLog('=== Initialization complete ===');
}

// Global error logger to surface uncaught errors
window.addEventListener('error', (e) => {
  debugLog('Uncaught error: ' + e.message);
});

// Boot: run init after DOM is ready
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', init);
} else {
  init();
}
window.addEventListener('load', init);