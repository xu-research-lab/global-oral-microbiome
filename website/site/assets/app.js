/* Global Oral Microbiome Atlas — front-end. Vanilla JS, no dependencies.
   Loads static JSON aggregates from ./data/ produced by precompute/build_data.py. */
'use strict';

const NON_WESTERN = new Set(['Asia', 'Africa', 'South America']);
const DATA = {};
const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));
const fmt = n => n == null ? '—' : n.toLocaleString('en-US');
const pct = x => (x == null ? '—' : (x * 100).toFixed(x >= 0.1 ? 1 : 2) + '%');
const cap = s => s ? s.charAt(0).toUpperCase() + s.slice(1) : s;

async function loadAll() {
  const names = ['overview', 'studies', 'taxa', 'refgap'];
  const results = await Promise.all(names.map(n =>
    fetch(`data/${n}.json`).then(r => {
      if (!r.ok) throw new Error(`${n}.json ${r.status}`);
      return r.json();
    })
  ));
  names.forEach((n, i) => DATA[n] = results[i]);
}

/* ---------- tab switching (with hash deep-links) ---------- */
function activateTab(name) {
  if (!$(`#tab-${name}`)) name = 'overview';
  $$('#tabs button').forEach(x => x.classList.toggle('active', x.dataset.tab === name));
  $$('section.tab').forEach(s => s.classList.toggle('active', s.id === `tab-${name}`));
}
function initTabs() {
  $('#tabs').addEventListener('click', e => {
    const b = e.target.closest('button[data-tab]');
    if (!b) return;
    location.hash = b.dataset.tab;
  });
  window.addEventListener('hashchange', () => activateTab(location.hash.slice(1)));
  activateTab(location.hash.slice(1) || 'overview');
}

/* ---------- bar chart helper ---------- */
function barChart(el, entries, { warnKeys = null, max = null, fmtVal = fmt } = {}) {
  const top = max ? entries.slice(0, max) : entries;
  const m = Math.max(...top.map(e => e[1]), 1);
  el.innerHTML = '<div class="bars">' + top.map(([label, val]) => {
    const w = (val / m * 100).toFixed(1);
    const warn = warnKeys && warnKeys.has(label) ? ' warn' : '';
    return `<div class="bar-row"><div class="lab" title="${label}">${label}</div>` +
      `<div class="bar-track"><div class="bar-fill${warn}" style="width:${w}%"></div></div>` +
      `<div class="val">${fmtVal(val)}</div></div>`;
  }).join('') + '</div>';
}

/* ---------- overview ---------- */
function renderOverview() {
  const o = DATA.overview, t = o.totals;
  const tiles = [
    ['samples', 'QC-passed samples'], ['studies', 'BioProjects'],
    ['subjects', 'Subjects'], ['genera', 'Genus-level taxa'],
    ['sites', 'Oral sites'], ['continents', 'Continents'],
  ];
  $('#stat-tiles').innerHTML = tiles.map(([k, l]) =>
    `<div class="stat"><div class="n">${fmt(t[k])}</div><div class="l">${l}</div></div>`).join('');
  $('#foot-stats').textContent = `${fmt(t.samples)} samples · ${fmt(t.studies)} studies · ${fmt(t.genera)} taxa`;

  barChart($('#bars-continent'), Object.entries(o.by_continent), { warnKeys: NON_WESTERN });
  barChart($('#bars-site'), Object.entries(o.by_site).map(([k, v]) => [cap(k), v]));
  barChart($('#bars-amplicon'), Object.entries(o.by_amplicon), { max: 8 });
  barChart($('#bars-disease'), Object.entries(o.disease_case_counts).filter(([k]) => k !== 'Other'), { max: 10 });
}

/* ---------- studies ---------- */
const studyState = { sort: 'n_samples', dir: -1 };
function studyColumns() {
  return [
    { k: 'bioproject', t: 'Accession', num: false },
    { k: 'title', t: 'Study', num: false },
    { k: 'continents', t: 'Geography', num: false, arr: true },
    { k: 'sites', t: 'Sites', num: false, arr: true },
    { k: 'amplicons', t: 'Amplicon', num: false, arr: true },
    { k: 'n_samples', t: 'Samples', num: true },
    { k: 'n_subjects', t: 'Subjects', num: true },
    { k: 'disease', t: 'Disease', num: false },
  ];
}
function initStudyFilters() {
  const s = DATA.studies.studies;
  const uniq = key => [...new Set(s.flatMap(r => r[key]))].filter(Boolean).sort();
  const fill = (sel, vals) => { $(sel).insertAdjacentHTML('beforeend',
    vals.map(v => `<option value="${v}">${cap(v)}</option>`).join('')); };
  fill('#study-continent', uniq('continents'));
  fill('#study-site', uniq('sites'));
  fill('#study-amplicon', uniq('amplicons'));
  ['#study-search', '#study-continent', '#study-site', '#study-amplicon', '#study-disease']
    .forEach(id => $(id).addEventListener('input', renderStudies));
}
function renderStudies() {
  const q = $('#study-search').value.toLowerCase().trim();
  const fc = $('#study-continent').value, fs = $('#study-site').value,
    fa = $('#study-amplicon').value, fd = $('#study-disease').value;
  let rows = DATA.studies.studies.filter(r =>
    (!q || r.bioproject.toLowerCase().includes(q) || (r.title || '').toLowerCase().includes(q)) &&
    (!fc || r.continents.includes(fc)) && (!fs || r.sites.includes(fs)) &&
    (!fa || r.amplicons.includes(fa)) && (fd !== '__disease' || r.is_disease_study));
  const sk = studyState.sort;
  rows.sort((a, b) => {
    let x = a[sk], y = b[sk];
    if (Array.isArray(x)) { x = x.length; y = y.length; }
    if (typeof x === 'string') return studyState.dir * x.localeCompare(y);
    return studyState.dir * ((x || 0) - (y || 0));
  });
  const cols = studyColumns();
  const head = '<thead><tr>' + cols.map(c => {
    const arr = studyState.sort === c.k ? `<span class="arrow">${studyState.dir < 0 ? '▼' : '▲'}</span>` : '';
    return `<th class="sortable" data-k="${c.k}">${c.t} ${arr}</th>`;
  }).join('') + '</tr></thead>';
  const body = '<tbody>' + rows.map(r => '<tr>' + cols.map(c => {
    if (c.k === 'bioproject') return `<td><a href="https://www.ncbi.nlm.nih.gov/bioproject/${r.bioproject}" target="_blank" rel="noopener">${r.bioproject}</a></td>`;
    if (c.k === 'title') return `<td>${r.title || '<span class="muted">—</span>'}</td>`;
    if (c.arr) return `<td>${r[c.k].map(v => `<span class="chip">${cap(v)}</span>`).join('') || '—'}</td>`;
    if (c.k === 'disease') return `<td>${r.is_disease_study
      ? r.diseases.filter(d => d !== 'Other').slice(0, 3).map(d => `<span class="chip warn">${d}</span>`).join('') || '<span class="chip warn">disease</span>'
      : '<span class="chip acc">healthy</span>'}</td>`;
    if (c.num) return `<td class="num">${fmt(r[c.k])}</td>`;
    return `<td>${r[c.k]}</td>`;
  }).join('') + '</tr>').join('') + '</tbody>';
  $('#study-table').innerHTML = head + body;
  $('#study-count').textContent = `${rows.length} of ${DATA.studies.studies.length} studies`;
  $$('#study-table th.sortable').forEach(th => th.addEventListener('click', () => {
    const k = th.dataset.k;
    if (studyState.sort === k) studyState.dir *= -1;
    else { studyState.sort = k; studyState.dir = (k === 'bioproject' || k === 'title') ? 1 : -1; }
    renderStudies();
  }));
}

/* ---------- taxa ---------- */
const taxaState = { sort: 'prev', dir: -1, limit: 250 };
function renderTaxaTable() {
  const q = $('#taxa-search').value.toLowerCase().trim();
  let rows = DATA.taxa.taxa;
  if (q) rows = rows.filter(t =>
    t.genus.toLowerCase().includes(q) || (t.family || '').toLowerCase().includes(q) ||
    (t.phylum || '').toLowerCase().includes(q));
  const total = rows.length;
  rows = rows.slice().sort((a, b) => taxaState.dir * ((a[taxaState.sort] || 0) - (b[taxaState.sort] || 0)));
  const shown = rows.slice(0, taxaState.limit);
  const head = '<thead><tr>' +
    '<th>Genus</th><th>Phylum</th>' +
    `<th class="sortable" data-k="prev">Prevalence ${taxaState.sort === 'prev' ? '▼' : ''}</th>` +
    `<th class="sortable" data-k="mean">Mean abund. ${taxaState.sort === 'mean' ? '▼' : ''}</th>` +
    '</tr></thead>';
  const body = '<tbody>' + shown.map(t =>
    `<tr data-id="${t.id}"><td><strong>${t.genus}</strong>${t.family ? `<br><span class="muted" style="font-size:11px">${t.family}</span>` : ''}</td>` +
    `<td><span class="chip">${t.phylum || '—'}</span></td>` +
    `<td class="num">${pct(t.prev)}</td><td class="num">${pct(t.mean)}</td></tr>`).join('') + '</tbody>';
  $('#taxa-table').innerHTML = head + body;
  $('#taxa-count').textContent = `${fmt(total)} taxa${q ? ' matched' : ''}` +
    (total > taxaState.limit ? ` · showing top ${taxaState.limit} by ${taxaState.sort === 'prev' ? 'prevalence' : 'abundance'}` : '');
  $$('#taxa-table tbody tr').forEach(tr => tr.addEventListener('click', () => showTaxon(+tr.dataset.id)));
  $$('#taxa-table th.sortable').forEach(th => th.addEventListener('click', () => {
    taxaState.sort = th.dataset.k; taxaState.dir = -1; renderTaxaTable();
  }));
}
function groupBars(el, groups, labels) {
  // groups: {idx: [mean, prev]}; render prevalence + mean bars per label.
  const rows = labels.map((lab, i) => [lab, groups[i]]).filter(([, g]) => g);
  if (!rows.length) { el.innerHTML = '<p class="muted">No data.</p>'; return; }
  const maxMean = Math.max(...rows.map(([, g]) => g[0]), 1e-9);
  el.innerHTML = '<div class="bars">' + rows.map(([lab, g]) => {
    const [mean, prev] = g;
    return `<div class="bar-row"><div class="lab" title="${lab}">${cap(lab)}</div>` +
      `<div class="bar-track" style="background:transparent">` +
        `<div class="bar-fill" style="width:${(prev * 100).toFixed(1)}%;opacity:.85"></div>` +
      `</div><div class="val">${pct(prev)}</div></div>` +
      `<div class="bar-row"><div class="lab"></div>` +
      `<div class="bar-track" style="background:transparent">` +
        `<div class="bar-fill" style="width:${(mean / maxMean * 100).toFixed(1)}%;background:var(--pos)"></div>` +
      `</div><div class="val">${pct(mean)}</div></div>`;
  }).join('') + '</div>';
}
function showTaxon(id) {
  const t = DATA.taxa.taxa.find(x => x.id === id);
  if (!t) return;
  const m = DATA.taxa.meta;
  $('#taxon-placeholder').style.display = 'none';
  $('#taxon-detail').classList.add('show');
  $('#td-name').textContent = t.genus;
  $('#td-tax').textContent = t.tax;
  const keyed = obj => { const o = {}; for (const k in obj) o[+k] = obj[k]; return o; };
  groupBars($('#td-site'), keyed(t.by_site), m.sites);
  groupBars($('#td-continent'), keyed(t.by_continent), m.continents);
  groupBars($('#td-status'), keyed(t.by_status), m.status);
}

/* ---------- reference gap ---------- */
function renderRefgap() {
  const r = DATA.refgap;
  $('#refgap-thesis').innerHTML =
    'After adjusting for amplicon region, oral reads from <strong>Asia</strong> map to the GG2 reference at ' +
    '~half the odds (OR 0.48) and from <strong>Africa</strong> at ~one-third the odds (OR 0.31) of North America, ' +
    'while <strong>Western Europe is indistinguishable</strong> from North America. The unmapped (reference-gap) ' +
    'fraction roughly doubles in Asian and triples in African studies — diversity the reference is missing, not low quality.';

  // Forest plot on a log-OR axis spanning 0.2 .. 1.3
  const lo = Math.log(0.2), hi = Math.log(1.35);
  const xpos = or => ((Math.log(or) - lo) / (hi - lo) * 100);
  const one = xpos(1).toFixed(1);
  $('#forest').innerHTML = r.regression.map(row => {
    const l = xpos(Math.exp(row.ci[0])), rr = xpos(Math.exp(row.ci[1])), p = xpos(row.or);
    return `<div class="forest-row"><div class="fl">${row.continent}</div>` +
      `<div class="forest-track"><div class="forest-axis" style="left:${one}%"></div>` +
      `<div class="forest-ci" style="left:${l}%;width:${(rr - l).toFixed(1)}%"></div>` +
      `<div class="forest-pt ${row.sig ? 'sig' : ''}" style="left:${p}%"></div></div>` +
      `<div class="fr">OR ${row.or}${row.sig ? ` · p=${row.p < 1e-3 ? row.p.toExponential(0) : row.p}` : ' · ns'}</div></div>`;
  }).join('') + `<div class="forest-row"><div class="fl muted" style="font-weight:400">reference</div>` +
    `<div class="forest-track"><div class="forest-axis" style="left:${one}%"></div>` +
    `<div class="forest-pt" style="left:${one}%;background:var(--muted)"></div></div>` +
    `<div class="fr">North America (OR 1)</div></div>`;

  const rows = r.summary;
  const head = '<thead><tr><th>Subset</th><th>Continent</th><th>Samples</th><th>Studies</th>' +
    '<th>Median mapping</th><th>Mean unmapped</th><th>% below 0.30 floor</th></tr></thead>';
  const body = '<tbody>' + rows.map(row => {
    const warn = NON_WESTERN.has(row.continent) ? ' style="color:var(--warn)"' : '';
    return `<tr><td>${row.subset}</td><td${warn}>${row.continent}</td>` +
      `<td class="num">${fmt(row.n_samples)}</td><td class="num">${fmt(row.n_studies)}</td>` +
      `<td class="num">${row.median_mapping}</td><td class="num">${row.mean_unmapped_pct}%</td>` +
      `<td class="num">${row.pct_below_floor}%</td></tr>`;
  }).join('') + '</tbody>';
  $('#refgap-table').innerHTML = head + body;
}

/* ---------- download tab ---------- */
const DL_FILES = [
  ['biom_oral_qc.biom', '53 MB', 'Oral feature table (OTU), QC-passed, non-deduplicated'],
  ['biom_oral_qc_dedup.biom', '35 MB', 'Oral feature table, one sample per subject × site'],
  ['biom_all_qc.biom', '70 MB', 'All body sites, QC-passed'],
  ['genus_table_relabund_all_qc.qza', '75 MB', 'Genus-level relative abundance (QIIME2 artifact)'],
  ['metadata_oral_qc.tsv', '31 MB', 'Harmonized oral metadata (matches biom_oral_qc)'],
  ['metadata_oral_qc_dedup.tsv', '19 MB', 'Harmonized oral metadata, deduplicated'],
  ['metadata_all_qc.tsv', '39 MB', 'Harmonized all-site metadata'],
];
const DICT = [
  ['sample_id', 'Unique sample identifier (matches BIOM columns)'],
  ['bioproject_id', 'Source BioProject / study accession'],
  ['subject_id', 'Subject identifier within a study'],
  ['site_cleaned', 'Harmonized oral site (saliva, tooth, tongue, oral mucosa, oral swab, gingiva)'],
  ['continent_group', 'Continent of sample origin'],
  ['country_cleaned', 'Harmonized country'],
  ['sequencing_region_cleaned', 'Harmonized 16S amplicon region (V4, V3-V4, …)'],
  ['age_cleaned_group', 'Age band (0-2, 3-17, 18-60, >60)'],
  ['sex_cleaned', 'Harmonized sex (F/M)'],
  ['status_cleaned', 'Health status — "HC" = healthy control, else disease name'],
  ['sample_cr_mapping_rate', 'Fraction of reads mapping to GG2 reference (99% id)'],
  ['Shannon', 'Shannon alpha diversity'],
];
function renderDownload() {
  $('#dl-list').innerHTML = DL_FILES.map(([f, s, d]) =>
    `<div class="dl-item"><div><code>${f}</code><br><span class="muted">${d}</span></div><div class="dl-size">${s}</div></div>`).join('');
  $('#dict-table').innerHTML = '<thead><tr><th>Field</th><th>Description</th></tr></thead><tbody>' +
    DICT.map(([f, d]) => `<tr><td><code>${f}</code></td><td>${d}</td></tr>`).join('') + '</tbody>';
}

/* ---------- boot ---------- */
(async function () {
  initTabs();
  try {
    await loadAll();
  } catch (e) {
    document.querySelector('main').innerHTML =
      `<div class="card"><h2>Could not load data</h2><p class="muted">${e.message}. ` +
      `Serve this folder over HTTP (e.g. <code>python3 -m http.server</code>) rather than opening the file directly.</p></div>`;
    return;
  }
  renderOverview();
  initStudyFilters(); renderStudies();
  renderTaxaTable();
  $('#taxa-search').addEventListener('input', renderTaxaTable);
  renderRefgap();
  renderDownload();
})();
