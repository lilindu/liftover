// Global state
let blocksAB = null;
let blocksBA = null;
const APP_VERSION = 'v1.02 (2025-12-25)';

const GENOME_PAIRS = {
  'pombase_leupold': {
    name: 'PomBase ↔ Leupold Consensus',
    genomeA: 'PomBase',
    genomeB: 'Leupold',
    fileAB: 'data/pombase_leupold/A_to_B.blocks.tsv',
    fileBA: 'data/pombase_leupold/B_to_A.blocks.tsv'
  },
  'pombase_dy47073': {
    name: 'PomBase ↔ DY47073',
    genomeA: 'PomBase',
    genomeB: 'DY47073',
    fileAB: 'data/pombase_dy47073/A_to_B.blocks.tsv',
    fileBA: 'data/pombase_dy47073/B_to_A.blocks.tsv'
  },
  'pombase_dy47071': {
    name: 'PomBase ↔ DY47071',
    genomeA: 'PomBase',
    genomeB: 'DY47071',
    fileAB: 'data/pombase_dy47071/A_to_B.blocks.tsv',
    fileBA: 'data/pombase_dy47071/B_to_A.blocks.tsv'
  }
};

let currentPairId = 'pombase_leupold';

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

function toZeroBased(posOneBased) {
  return posOneBased - 1;
}

function toOneBased(posZeroBased) {
  return posZeroBased + 1;
}

// Fetch default A→B blocks from server with cache-busting
async function loadDefaultAB() {
  const conf = GENOME_PAIRS[currentPairId];
  const statusAB = document.getElementById('statusAB');
  try {
    debugLog(`Fetching ${conf.genomeA}→${conf.genomeB} blocks...`);
    const res = await fetch(conf.fileAB + '?v=' + Date.now());
    if (!res.ok) throw new Error('HTTP ' + res.status);
    const txt = await res.text();
    blocksAB = parseTSV(txt);
    const count = Object.keys(blocksAB).reduce((s, c) => s + blocksAB[c].length, 0);
    if (statusAB) statusAB.textContent = `${conf.genomeA}→${conf.genomeB}: blocks file loaded (${count} blocks)`;
    debugLog(`${conf.genomeA}→${conf.genomeB} blocks file loaded: ${count} blocks`);
  } catch (e) {
    if (statusAB) statusAB.textContent = `${conf.genomeA}→${conf.genomeB}: failed (${e.message})`;
    debugLog(`${conf.genomeA}→${conf.genomeB} failed: ${e.message}`);
  }
}

// Fetch default B→A blocks from server with cache-busting
async function loadDefaultBA() {
  const conf = GENOME_PAIRS[currentPairId];
  const statusBA = document.getElementById('statusBA');
  try {
    debugLog(`Fetching ${conf.genomeB}→${conf.genomeA} blocks...`);
    const res = await fetch(conf.fileBA + '?v=' + Date.now());
    if (!res.ok) throw new Error('HTTP ' + res.status);
    const txt = await res.text();
    blocksBA = parseTSV(txt);
    const count = Object.keys(blocksBA).reduce((s, c) => s + blocksBA[c].length, 0);
    if (statusBA) statusBA.textContent = `${conf.genomeB}→${conf.genomeA}: blocks file loaded (${count} blocks)`;
    debugLog(`${conf.genomeB}→${conf.genomeA} blocks file loaded: ${count} blocks`);
  } catch (e) {
    if (statusBA) statusBA.textContent = `${conf.genomeB}→${conf.genomeA}: failed (${e.message})`;
    debugLog(`${conf.genomeB}→${conf.genomeA} failed: ${e.message}`);
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
    const posOneBased = parseInt(coords);
    if (!Number.isFinite(posOneBased) || posOneBased < 1) return null;
    const posZeroBased = toZeroBased(posOneBased);
    return { contig, start: posZeroBased, end: posZeroBased, isInterval: false };
  }
  const m = coords.match(/^(\d+)-(\d+)$/);
  if (!m) return null;
  const aOneBased = parseInt(m[1]);
  const bOneBased = parseInt(m[2]);
  if (!Number.isFinite(aOneBased) || !Number.isFinite(bOneBased)) return null;
  if (aOneBased < 1 || bOneBased < 1) return null;
  const startOneBased = Math.min(aOneBased, bOneBased);
  const endOneBased = Math.max(aOneBased, bOneBased);
  const startZeroBased = toZeroBased(startOneBased);
  const endZeroBased = toZeroBased(endOneBased);
  return { contig, start: startZeroBased, end: endZeroBased, isInterval: true };
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
      gaps.push({ type: 'GAP_A', sizeA: gapA, aStart: toOneBased(cur.endA), aEnd: toOneBased(next.startA - 1) });
    }
    const bEnd = mapPoint(cur, cur.endA - 1)[1];
    const bStartNext = mapPoint(next, next.startA)[1];
    if (targetStrand === '+') {
      const gapB = Math.max(0, bStartNext - bEnd - 1);
      const ovlB = Math.max(0, bEnd - bStartNext + 1);
      if (gapB > 0) gaps.push({ type: 'GAP_B', sizeB: gapB, bStart: toOneBased(bEnd + 1), bEnd: toOneBased(bStartNext - 1) });
      if (ovlB > 0) gaps.push({ type: 'OVERLAP_B', sizeB: ovlB, bStart: toOneBased(bStartNext), bEnd: toOneBased(bEnd) });
    } else {
      const gapB = Math.max(0, bEnd - bStartNext - 1);
      const ovlB = Math.max(0, bStartNext - bEnd + 1);
      if (gapB > 0) {
        const lo = bStartNext + 1;
        const hi = bEnd - 1;
        gaps.push({ type: 'GAP_B', sizeB: gapB, bStart: toOneBased(Math.min(lo, hi)), bEnd: toOneBased(Math.max(lo, hi)) });
      }
      if (ovlB > 0) {
        gaps.push({ type: 'OVERLAP_B', sizeB: ovlB, bStart: toOneBased(Math.min(bStartNext, bEnd)), bEnd: toOneBased(Math.max(bStartNext, bEnd)) });
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
  updateDirectionNote && updateDirectionNote('AB');
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
    if (!blocks) { out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t\t\t\t\tNO_CONTIG`); continue; }
    const attempt = stitchInterval(blocks, start, end);
    if (!attempt.ok) {
      out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t\t\t\t\tSTITCH_FAILED_${attempt.reason}\t`);
      continue;
    }
    const gaps = (findBlock(blocks, start)===findBlock(blocks, end)) ? [] : collectGaps(blocks, findBlock(blocks, start), findBlock(blocks, end));
    const status = (gaps.length === 0) ? ((findBlock(blocks, start)===findBlock(blocks, end))?'OK':'STITCHED_OK') : 'STITCHED_WITH_GAPS';
    out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t${attempt.contigB}\t${toOneBased(attempt.startB)}\t${toOneBased(attempt.endB)}\t${attempt.strand}\t${status}\t${formatGaps(gaps)}`);
  }
  document.getElementById('out').textContent = out.join('\n');
  renderTSVToTable(document.getElementById('out').textContent);
}

// Run liftover Leupold→PomBase for input lines (supports points and intervals)
function runLiftoverBA() {
  updateDirectionNote && updateDirectionNote('BA');
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
    if (!blocks) { out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t\t\t\t\tNO_CONTIG`); continue; }
    const attempt = stitchInterval(blocks, start, end);
    if (!attempt.ok) {
      out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t\t\t\t\tSTITCH_FAILED_${attempt.reason}\t`);
      continue;
    }
    const gaps = (findBlock(blocks, start)===findBlock(blocks, end)) ? [] : collectGaps(blocks, findBlock(blocks, start), findBlock(blocks, end));
    const status = (gaps.length === 0) ? ((findBlock(blocks, start)===findBlock(blocks, end))?'OK':'STITCHED_OK') : 'STITCHED_WITH_GAPS';
    out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t${attempt.contigB}\t${toOneBased(attempt.startB)}\t${toOneBased(attempt.endB)}\t${attempt.strand}\t${status}\t${formatGaps(gaps)}`);
  }
  document.getElementById('out').textContent = out.join('\n');
  renderTSVToTable(document.getElementById('out').textContent);
}

// Run round-trip PomBase→Leupold→PomBase for input lines (supports points and intervals)
function runRoundtrip() {
  updateDirectionNote && updateDirectionNote('RT');
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
    if (!blocks1) { out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t\t\t\t\t\t\tNO_CONTIG`); continue; }
    const attempt1 = stitchInterval(blocks1, start, end);
    if (!attempt1.ok) { out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t\t\t\t\t\t\tSTITCH_FAILED_${attempt1.reason}\t\t`); continue; }
    const gapsAB = (findBlock(blocks1, start)===findBlock(blocks1, end)) ? [] : collectGaps(blocks1, findBlock(blocks1, start), findBlock(blocks1, end));

    // Stage 2: Leupold→PomBase stitch
    const blocks2 = blocksBA[attempt1.contigB];
    if (!blocks2) { out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t${attempt1.contigB}\t${toOneBased(attempt1.startB)}\t${toOneBased(attempt1.endB)}\t\t\t\tNO_CONTIG_BA`); continue; }
    const attempt2 = stitchInterval(blocks2, attempt1.startB, attempt1.endB);
    if (!attempt2.ok) { out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t${attempt1.contigB}\t${toOneBased(attempt1.startB)}\t${toOneBased(attempt1.endB)}\t\t\t\tSTITCH_FAILED_BA_${attempt2.reason}\t${formatGaps(gapsAB)}\t`); continue; }
    const gapsBA = (findBlock(blocks2, attempt1.startB)===findBlock(blocks2, attempt1.endB)) ? [] : collectGaps(blocks2, findBlock(blocks2, attempt1.startB), findBlock(blocks2, attempt1.endB));

    const status = (attempt2.contigB === contig && Math.min(attempt2.startB, attempt2.endB) === start && Math.max(attempt2.startB, attempt2.endB) === end) ? 'PASS' : 'FAIL';
    out.push(`${contig}\t${toOneBased(start)}\t${toOneBased(end)}\t${attempt1.contigB}\t${toOneBased(attempt1.startB)}\t${toOneBased(attempt1.endB)}\t${attempt2.contigB}\t${toOneBased(Math.min(attempt2.startB, attempt2.endB))}\t${toOneBased(Math.max(attempt2.startB, attempt2.endB))}\t${status}\t${formatGaps(gapsAB)}\t${formatGaps(gapsBA)}`);
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

// Populate pair selector and handle changes
function initPairSelector() {
  const sel = document.getElementById('pairSelector');
  if (!sel) return;
  sel.innerHTML = '';
  for (const pid in GENOME_PAIRS) {
    const opt = document.createElement('option');
    opt.value = pid;
    opt.textContent = GENOME_PAIRS[pid].name;
    if (pid === currentPairId) opt.selected = true;
    sel.appendChild(opt);
  }
  sel.addEventListener('change', (e) => {
    switchPair(e.target.value);
  });
}

// Switch current pair and reload
async function switchPair(newId) {
  if (!GENOME_PAIRS[newId]) return;
  currentPairId = newId;
  debugLog('Switching pair to: ' + newId);
  blocksAB = null;
  blocksBA = null;
  
  // Update UI text immediately
  updateUIText();
  
  // Clear status
  const sAB = document.getElementById('statusAB');
  const sBA = document.getElementById('statusBA');
  if (sAB) sAB.textContent = 'A→B: loading...';
  if (sBA) sBA.textContent = 'B→A: loading...';
  
  // Clear results
  document.getElementById('out').textContent = '';
  const head = document.getElementById('outHead');
  const body = document.getElementById('outBody');
  if (head) head.innerHTML = '';
  if (body) body.innerHTML = '';
  
  // Load new blocks
  await loadDefaultAB();
  await loadDefaultBA();
}

// Update buttons and labels based on current pair
function updateUIText() {
  const conf = GENOME_PAIRS[currentPairId];
  const btnAB = document.getElementById('runAB');
  const btnBA = document.getElementById('runBA');
  const btnRT = document.getElementById('runRoundtrip');
  
  if (btnAB) btnAB.textContent = `Liftover ${conf.genomeA}→${conf.genomeB}`;
  if (btnBA) btnBA.textContent = `Liftover ${conf.genomeB}→${conf.genomeA}`;
  if (btnRT) btnRT.textContent = `Round-Trip ${conf.genomeA}→${conf.genomeB}→${conf.genomeA}`;
}

// Initialize: attach handlers and autoload blocks
async function init() {
  debugLog('=== Starting autoload initialization ===');
  debugLog('DOM readyState: ' + document.readyState);
  const vt = document.getElementById('versionText');
  if (vt) vt.textContent = APP_VERSION;

  initPairSelector();
  updateUIText();

  const sAB = document.getElementById('statusAB');
  const sBA = document.getElementById('statusBA');
  const conf = GENOME_PAIRS[currentPairId];
  
  if (sAB) sAB.textContent = `${conf.genomeA}→${conf.genomeB}: autoload starting`;
  if (sBA) sBA.textContent = `${conf.genomeB}→${conf.genomeA}: autoload starting`;

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
// Write version badge if present
const vb = document.getElementById('versionBadge');
if (vb) vb.textContent = APP_VERSION;

function updateDirectionNote(mode) {
  const conf = GENOME_PAIRS[currentPairId];
  const dn = document.getElementById('directionNote');
  if (!dn) return;
  if (mode === 'AB') dn.textContent = `Current direction: A = ${conf.genomeA} → B = ${conf.genomeB}`;
  else if (mode === 'BA') dn.textContent = `Current direction: A = ${conf.genomeB} → B = ${conf.genomeA}`;
  else if (mode === 'RT') dn.textContent = `Current direction: A = ${conf.genomeA} → B = ${conf.genomeB} → A = ${conf.genomeA}`;
}
