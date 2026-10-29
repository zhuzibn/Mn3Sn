/* Display helpers consume the MATLAB snapshot, never construct physical geometry. */
(function (root) {
  'use strict';
  const data = root.MN3SN_DATA;
  function snapshot(W, L) {
    const s = data.snapshots.find(s => s.W === W && s.L === L);
    if (!s) throw new Error('Unsupported static snapshot');
    return s;
  }
  function rc(s, i) { return [i % s.W + 1, Math.floor(i / s.W) + 1]; }
  function index(s, r, c) { return r - 1 + (c - 1) * s.W; }
  function position(r, c) { return [c - 1 + (r % 2 ? 0.5 : 0), (r - 1) * Math.sqrt(3) / 2]; }
  function group(s, i) {
    if (s.occupancy[i] !== 1) return 'empty';
    if (s.axes[i][1] < -0.5) return 'C';
    return s.axes[i][2] > 0 ? 'A' : 'B';
  }
  function reads(s, i, term, bc) {
    return data.offsets.map(([dr, dc], q) => {
      const [r, c] = rc(s, i), rr = r + dr, cc = c + dc;
      const outside = rr < 1 || rr > s.W || cc < 1 || cc > s.L;
      const wr = ((rr - 1) % s.W + s.W) % s.W + 1;
      const wc = ((cc - 1) % s.L + s.L) % s.L + 1;
      const j = outside && bc === 1 ? null : index(s, wr, wc);
      const raw = (term === 'exchange' ? s.exchange : s.dmi)[i][q];
      const weight = raw * s.scale[i];
      let status = 'active';
      if (s.occupancy[i] !== 1) status = 'empty receiver';
      else if (raw === 0) status = 'coefficient cut';
      else if (j === null) status = 'open edge';
      else if (s.occupancy[j] !== 1) status = 'empty source';
      return {q, name: data.shiftNames[q], dr, dc, rr, cc, j, raw, weight,
        outside, wrapped: outside && bc === 0, wrap: [Math.floor((rr - 1) / s.W), Math.floor((cc - 1) / s.L)],
        status, active: status === 'active'};
    });
  }
  function reverse(s, i, read, term, bc) {
    if (read.j === null) return 0;
    return reads(s, read.j, term, bc).filter(r => r.active && r.j === i).reduce((sum, r) => sum + r.weight, 0);
  }
  function exchangePairs(s, bc) {
    const pairs = new Map();
    for (let i = 0; i < s.occupancy.length; i++) {
      for (const read of reads(s, i, 'exchange', bc).filter(r => r.active)) {
        const a = Math.min(i, read.j), b = Math.max(i, read.j);
        if (a === b) continue;
        const key = a + ':' + b;
        if (!pairs.has(key)) pairs.set(key, {a, b, intoA: 0, intoB: 0});
        const bond = pairs.get(key);
        bond[i === a ? 'intoA' : 'intoB'] += read.weight;
      }
    }
    return [...pairs.values()].sort((x, y) => x.a - y.a || x.b - y.b);
  }
  root.Mn3Sn = {data, snapshot, rc, index, position, group, reads, reverse, exchangePairs};
})(typeof window !== 'undefined' ? window : globalThis);
