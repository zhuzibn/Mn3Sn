/* Display only: consume the verified production bond table, including image offsets. */
(() => {
  'use strict';
  const $ = id => document.getElementById(id), data = window.PERIODIC_DATA;
  const colors = ['#0072b2', '#d55e00', '#009e73'], names = ['a', 'b', 'c'];
  let selected = 0, neighbor = null, pan = [0, 0], drag = null;
  const graph = () => data.snapshots.find(s => s.tag === $('size').value);
  $('size').innerHTML = data.snapshots.map(s => `<option value="${s.tag}">${s.W} × ${s.L} / ${s.W / 4} × ${s.L / 2} motif repeats</option>`).join('');
  $('size').value = '8x8';
  const position = (r, c) => [c - 1 + (r % 2 ? .5 : 0), (r - 1) * Math.sqrt(3) / 2];
  function resetSites() {
    const s = graph(); selected = 0; neighbor = null; pan = [0, 0]; $('zoom').value = 100;
    $('receiver').innerHTML = s.sites.map((v, i) => `<option value="${i}">${i}: ${names[v.cls]} / (${v.row}, ${v.col})</option>`).join('');
    render();
  }
  function render() {
    const s = graph(), site = s.sites[selected], reads = s.reads.filter(b => b.i === selected);
    const scale = Math.min(650 / (s.L + 2), 410 / (s.W * Math.sqrt(3) / 2 + 2));
    const point = (r, c) => { const p = position(r, c); return [90 + scale * (p[0] + 1), 50 + scale * (p[1] + 1)]; };
    const xy = i => point(s.sites[i].row, s.sites[i].col);
    const term = $('term').value;
    const seamCount = s.reads.filter(b => b.image_row || b.image_col).length / 2;
    $('stats').innerHTML = [[s.sites.length, 'unique spins'], [s.sites.length / 3, 'spins per class'], [4, 'neighbors per spin'], [seamCount, 'unique seam bonds']].map(([v, label]) => `<div class="stat"><strong>${v}</strong><span>${label}</span></div>`).join('');
    $('heading').textContent = $('view').value === 'cell' ? 'Complete repeating cells' : `${term === 'K' ? 'Onsite term' : 'Sources'} for spin ${selected}`;
    $('details').innerHTML = `<dl><dt>Occupied ID / class</dt><dd>${selected} / ${names[site.cls]}</dd><dt>Array row, column</dt><dd>(${site.row}, ${site.col})</dd><dt>4 × 2 motif translation</dt><dd>(${Math.floor((site.row - 1) / 4)}, ${Math.floor((site.col - 1) / 2)})</dd><dt>Array linear index (MATLAB)</dt><dd>${site.row + (site.col - 1) * s.W}</dd><dt>Neighbors by class</dt><dd>${names.map((name, c) => `${name}: ${reads.filter(b => s.sites[b.j].cls === c).length}`).join(' · ')}</dd></dl>`;
    $('formula').textContent = term === 'A' ? 'Bex,i = −(J/μ) Σ Aij mj; J/μ = 100.812 T' : term === 'C' ? 'BDMI,i = (D/μ) Σ Cij (0,mjz,−mjy); D = 0.833 meV' : 'Bani,i = (2K/μ)(ei·mi)ei; K = 0.196 meV at every site';
    let markup = `<g transform="translate(${pan[0]} ${pan[1]}) translate(430 270) scale(${$('zoom').value / 100}) translate(-430 -270)">`;
    const tl = point(.4, .4), br = point(s.W + .6, s.L + .6);
    markup += `<rect x="${tl[0]}" y="${tl[1]}" width="${br[0] - tl[0]}" height="${br[1] - tl[1]}" fill="none" stroke="#879b9c" stroke-width="1.5"/>`;
    // Motif frames are index bookkeeping, not crystallographic unit cells.
    if ($('view').value === 'cell') {
      for (let r = 1; r <= s.W; r += 4) for (let c = 1; c <= s.L; c += 2) {
        const p = point(r, c);
        markup += `<rect x="${p[0] - scale * .65}" y="${p[1] - scale * .4}" width="${scale * 2}" height="${scale * 4 * Math.sqrt(3) / 2}" fill="none" stroke="#ccd1ca" stroke-dasharray="3 5"/>`;
      }
    }
    if (term !== 'K') for (const b of s.reads) {
      const wrapped = Boolean(b.image_row || b.image_col), active = b.i === selected;
      if (!active && ($('view').value !== 'cell' || (wrapped && !$('seams').checked))) continue;
      if (!active && !wrapped && b.i > b.j) continue;
      const source = s.sites[b.i], p = xy(b.i), q = point(source.row + b.dr, source.col + b.dc);
      const stroke = active ? (term === 'C' && b.C < 0 ? '#b04642' : '#226a9b') : '#c8ceca';
      markup += `<line x1="${p[0]}" y1="${p[1]}" x2="${q[0]}" y2="${q[1]}" stroke="${stroke}" stroke-width="${active ? (b.j === neighbor ? 4 : 2.5) : 1}" ${wrapped ? 'stroke-dasharray="5 4"' : ''}/>`;
      if (wrapped && active) markup += `<circle cx="${q[0]}" cy="${q[1]}" r="7" fill="#fffdf8" stroke="${colors[s.sites[b.j].cls]}" stroke-width="2"/><text x="${q[0] + 10}" y="${q[1] + 4}" font-size="11">${b.j} image</text>`;
    }
    for (let i = 0; i < s.sites.length; i++) {
      const p = xy(i), v = s.sites[i];
      markup += `<g class="site" data-site="${i}" role="button" tabindex="0" aria-label="Spin ${i}, class ${names[v.cls]}, row ${v.row}, column ${v.col}"><circle cx="${p[0]}" cy="${p[1]}" r="${i === selected ? 10 : 6}" fill="${colors[v.cls]}" stroke="${i === selected ? '#b76b10' : i === neighbor ? '#18313b' : '#fffdf8'}" stroke-width="${i === selected || i === neighbor ? 3 : 1}"/>${$('labels').checked ? `<text x="${p[0] + 9}" y="${p[1] - 9}" font-size="10" fill="#18313b">${i}</text>` : ''}</g>`;
    }
    markup += `<text x="${tl[0]}" y="${br[1] + 25}" font-size="12">columns → / rows ↓ · schematic display units</text></g>`;
    $('graph').innerHTML = markup;
    $('bonds').innerHTML = term === 'K' ? '<tr><td colspan="7">Onsite anisotropy has no neighbor read. All occupied sites retain K = 0.196 meV.</td></tr>' : reads.map(b => {
      const other = s.sites[b.j], reverse = s.reads.find(v => v.i === b.j && v.j === selected);
      return `<tr class="${neighbor === b.j ? 'selected' : ''}" data-neighbor="${b.j}"><td><button data-neighbor="${b.j}">${b.j}: ${names[other.cls]}</button></td><td>(${other.row}, ${other.col})</td><td>(${b.dr}, ${b.dc})</td><td>(${b.image_row}, ${b.image_col})</td><td>${b.A}</td><td>${b.C > 0 ? '+' : ''}${b.C}</td><td>${reverse.A} / ${reverse.C}</td></tr>`;
    }).join('');
    $('status').textContent = `${s.tag}: ${s.sites.length} sites, ${s.reads.length} directed bonds. Data: MATLAB impulse export plus independent stencil check.`;
    $('graph').dataset.selected = selected;
  }
  $('size').onchange = resetSites;
  $('receiver').onchange = () => { selected = Number($('receiver').value); neighbor = null; render(); };
  for (const id of ['view', 'term', 'labels', 'seams', 'zoom']) $(id).oninput = render;
  $('reset').onclick = () => { pan = [0, 0]; $('zoom').value = 100; render(); };
  function choose(event) {
    const node = event.target.closest('[data-site]');
    if (!node) return;
    selected = Number(node.dataset.site); neighbor = null; $('receiver').value = selected; render();
  }
  $('graph').onclick = choose;
  $('graph').onkeydown = e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); choose(e); } };
  $('bonds').onclick = e => { const node = e.target.closest('[data-neighbor]'); if (node) { neighbor = Number(node.dataset.neighbor); render(); } };
  $('graph').onpointerdown = e => { if (!e.target.closest('[data-site]')) { drag = [e.clientX, e.clientY, ...pan]; $('graph').setPointerCapture(e.pointerId); } };
  $('graph').onpointermove = e => { if (drag) { const ratio = 860 / $('graph').getBoundingClientRect().width; pan = [drag[2] + (e.clientX - drag[0]) * ratio, drag[3] + (e.clientY - drag[1]) * ratio]; render(); } };
  $('graph').onpointerup = $('graph').onpointercancel = () => { drag = null; };
  resetSites();
})();
