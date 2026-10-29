(function () {
  'use strict';
  const M = window.Mn3Sn, data = M.data, $ = id => document.getElementById(id);
  const isCell = document.body.dataset.view === 'cell';
  const colors = {A:'#226a9b', B:'#b04642', C:'#38796b', empty:'#b3b9b4'};
  const params = new URLSearchParams(location.search);
  const requested = (params.get('size') || '5x6').split('x').map(Number);
  let s = data.snapshots.find(x => x.W === requested[0] && x.L === requested[1]) || M.snapshot(5,6);
  let selected = Math.max(0,Math.min(s.W*s.L-1,Math.trunc(Number(params.get('site'))) || 0));
  let readQ = -1, pan = [0,0], drag = null, moved = false;
  const esc = value => String(value).replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
  const pair = i => i === null ? 'outside' : '('+M.rc(s,i).join(', ')+')';
  const vec = v => '('+v.map(x => Math.abs(x)<1e-14?'0':x.toFixed(4)).join(', ')+')';
  const sci = x => x === 0 ? '0' : x.toExponential(4);
  const signed = x => x > 0 ? '+'+x : String(x);
  const term = () => isCell ? 'exchange' : $('term').value;
  const bc = () => isCell ? 1 : Number($('boundary').value);
  const activeReads = () => term()==='anisotropy' ? [] : M.reads(s,selected,term(),bc());
  function query() { return '?size='+s.W+'x'+s.L+'&site='+selected; }
  function setSize(W,L) {
    s=M.snapshot(W,L); selected=Math.min(selected,s.W*s.L-1); readQ=-1; pan=[0,0]; $('zoom').value=100;
    syncControls(); render();
  }
  function syncControls() {
    const value=s.W+'x'+s.L;
    const old=$('preset').querySelector('[data-custom]'); if(old) old.remove();
    if (![...$('preset').options].some(o=>o.value===value)) {
      const option=new Option(value.replace('x',' × ')+' / repeated patch',value); option.dataset.custom='true'; $('preset').add(option);
    }
    $('preset').value=value;
    if(isCell) {
      $('repeat-w').value=s.W;
      const old=$('repeat-l').querySelector('[data-custom]'); if(old) old.remove();
      if(s.L===5) {const option=new Option('5 columns / legacy edge',5);option.dataset.custom='true';$('repeat-l').add(option);}
      $('repeat-l').value=s.L;
    }
    $('target').innerHTML=s.occupancy.map((v,i)=>'<option value="'+i+'">'+pair(i)+' · '+(v===1?'class '+M.group(s,i):'empty')+'</option>').join('');
    $('target').value=selected;
  }
  function selectSite(i) {selected=i;readQ=-1;$('target').value=i;render();}
  function render() {
    const reads=activeReads(), occupied=s.occupancy.filter(x=>x===1).length;
    const probe=s.probes.find(p=>p.bc===bc());
    const values=isCell ? [[occupied,'occupied Mn entries'],[s.W*s.L-occupied,'empty slots'],[probe.exchangePairs,'exchange connections']]
      : [[reads.filter(r=>r.active).length,term()==='anisotropy'?'neighbor reads / onsite term':'active neighbor reads'],[probe.unequalPairs,'unequal exchange pairs'],[occupied,'occupied Mn entries']];
    $('stats').innerHTML=values.map(([v,l])=>'<div class="stat"><strong>'+v+'</strong><span>'+l+'</span></div>').join('');
    $('legend').innerHTML=['A','B','C'].map((g,i)=>'<span><i class="swatch" style="background:'+colors[g]+'"></i>Class '+g+' / θ = '+[30,150,270][i]+'°</span>').join('')+'<span>○ empty slot</span><span>Gold ring: selected entry</span><span>θ measured from spin-space +z toward +y</span>';
    $('other-view').href=(isCell?'structure-exchange.html':'unit-cell.html')+query();
    if(isCell) $('inspect-link').href='structure-exchange.html'+query();
    details(); renderSVG(); renderArray(); sourceExcerpts();
    if(!isCell) {renderTable(); renderFormula(); renderUnequalPairs();}
    $('status').textContent=s.W+' × '+s.L+' static snapshot · '+(bc()?'open reads':'legacy wrapped reads')+' · '+(s.occupancy[selected]===1?'selected Mn ':'selected empty slot ')+pair(selected);
  }
  function renderUnequalPairs() {
    $('unequal-exchange').hidden=term()!=='exchange';
    if(term()!=='exchange')return;
    const all=M.exchangePairs(s,bc()), unequal=all.filter(p=>Math.abs(p.intoA-p.intoB)>data.tolerance);
    const read=activeReads().find(r=>r.q===readQ);
    const selectedKey=read?.active?[Math.min(selected,read.j),Math.max(selected,read.j)].join(':'):'';
    const key=p=>p.a+':'+p.b, chosen=unequal.find(p=>key(p)===selectedKey);
    const weight=x=>x===0.5?'½':String(x);
    const shade=h=>h===1?'#a13f37':h===0.5?'#226a9b':'#38796b';
    $('unequal-count').textContent=unequal.length+' unequal / '+all.length+' total pairs';
    $('unequal-selection').textContent=chosen?'Selected pair '+(unequal.indexOf(chosen)+1)+': '+pair(chosen.a)+' ↔ '+pair(chosen.b)+'. The site inspector below follows this selection.':unequal.length?'Every unequal pair is shown. Select a number or table button to inspect it.':'No unequal exchange pairs in this snapshot.';
    const scale=Math.min(510/(s.L-.5),335/((s.W-1)*Math.sqrt(3)/2));
    const point=i=>{const p=M.position(...M.rc(s,i));return [310-(s.L-.5)*scale/2+p[0]*scale,250-(s.W-1)*Math.sqrt(3)/2*scale/2+p[1]*scale];};
    const font=Math.max(9,Math.min(16,scale*.18)),radius=Math.min(13,scale*.14);
    let svg='<title>Unequal directed exchange weights, '+s.W+' by '+s.L+'</title><desc>Numbered orange links correspond to the directional weights in the adjacent table. Node numbers show the receiver scaling factor.</desc>';
    svg+='<text x="20" y="28" font-size="14" fill="#53696f">columns → · rows ↓ · '+(bc()===1?'open patch':'wrapped endpoints; legacy cuts retained')+'</text>';
    for(const p of all) {
      const a=point(p.a),b=point(p.b);
      svg+='<line x1="'+a[0]+'" y1="'+a[1]+'" x2="'+b[0]+'" y2="'+b[1]+'" stroke="#c8d0ca" stroke-width="2"/>';
    }
    unequal.forEach((p,n)=>{
      const a=point(p.a),b=point(p.b),x=(a[0]+b[0])/2,y=(a[1]+b[1])/2,active=key(p)===selectedKey;
      const label='Pair '+(n+1)+': '+pair(p.b)+' to '+pair(p.a)+' weight '+p.intoA+' J; reverse '+p.intoB+' J';
      svg+='<g class="unequal-bond" data-exchange-pair="'+key(p)+'" role="button" tabindex="0" aria-pressed="'+active+'" aria-label="'+label+'"><title>'+label+'</title><line x1="'+a[0]+'" y1="'+a[1]+'" x2="'+b[0]+'" y2="'+b[1]+'" stroke="transparent" stroke-width="20"/><line x1="'+a[0]+'" y1="'+a[1]+'" x2="'+b[0]+'" y2="'+b[1]+'" stroke="'+(active?'#18313b':'#b76b10')+'" stroke-width="'+(active?6:3.5)+'"/><circle cx="'+x+'" cy="'+y+'" r="'+(font*.8)+'" fill="'+(active?'#18313b':'#fff4dc')+'" stroke="#b76b10"/><text x="'+x+'" y="'+(y+font*.35)+'" text-anchor="middle" font-size="'+font+'" fill="'+(active?'white':'#704109')+'">'+(n+1)+'</text></g>';
    });
    for(let i=0;i<s.occupancy.length;i++) {
      const [x,y]=point(i),empty=s.occupancy[i]!==1,h=s.scale[i];
      const active=chosen&&(i===chosen.a||i===chosen.b);
      svg+='<g><title>'+pair(i)+(empty?' empty':' receiver h = '+h)+'</title>';
      if(active)svg+='<circle cx="'+x+'" cy="'+y+'" r="'+(radius+5)+'" fill="none" stroke="#b76b10" stroke-width="2"/>';
      svg+='<circle cx="'+x+'" cy="'+y+'" r="'+(empty?radius*.55:radius)+'" fill="'+(empty?'#fffdf8':shade(h))+'" stroke="'+(empty?'#abb6af':shade(h))+'" '+(empty?'stroke-dasharray="3 3"':'')+'/>';
      if(!empty)svg+='<text x="'+x+'" y="'+(y+font*.35)+'" text-anchor="middle" font-size="'+font+'" fill="white">'+weight(h)+'</text>';
      svg+='<text x="'+x+'" y="'+(y+radius+font+3)+'" text-anchor="middle" font-size="'+font+'" fill="'+(empty?'#7d8a84':'#18313b')+'" stroke="#fffdf8" stroke-width="4" paint-order="stroke">'+M.rc(s,i).join(',')+'</text></g>';
    }
    $('unequal-map').setAttribute('viewBox','0 0 620 500');
    $('unequal-map').innerHTML=svg;
    const scales=[...new Set(s.scale.filter((h,i)=>s.occupancy[i]===1))].sort((a,b)=>b-a);
    $('unequal-legend').innerHTML=scales.map(h=>'<span><i class="swatch" style="background:'+shade(h)+'"></i>Receiver h = '+weight(h)+'</span>').join('')+'<span>○ empty slot</span>';
    $('unequal-table').querySelector('tbody').innerHTML=unequal.map((p,n)=>'<tr class="'+(key(p)===selectedKey?'selected':'')+'"><td><button type="button" data-exchange-pair="'+key(p)+'" aria-pressed="'+(key(p)===selectedKey)+'" aria-label="Inspect unequal pair '+(n+1)+'">'+(n+1)+'</button></td><td class="pair-directions">'+pair(p.b)+' → '+pair(p.a)+'<br>'+pair(p.a)+' → '+pair(p.b)+'</td><td>'+weight(p.intoA)+'<br>'+weight(p.intoB)+'</td></tr>').join('');
  }
  function selectExchangePair(value) {
    const [a,b]=value.split(':').map(Number);
    let receiver=a,read=M.reads(s,a,'exchange',bc()).find(r=>r.active&&r.j===b);
    if(!read){receiver=b;read=M.reads(s,b,'exchange',bc()).find(r=>r.active&&r.j===a);}
    if(!read)return;
    selected=receiver;readQ=read.q;$('target').value=selected;render();
  }
  function details() {
    const [r,c]=M.rc(s,selected), g=M.group(s,selected), xy=M.position(r,c);
    $('site-heading').textContent=(g==='empty'?'Empty entry ':'Mn entry ')+pair(selected);
    const entries=[['Array identity','row '+r+', column '+c+' · MATLAB linear index '+(selected+1)],
      ['Code class',g==='empty'?'No spin / atomtype = 3':g+' / atomtype = 1'],
      ['Illustration coordinates',vec(xy)+' · arbitrary units'],
      ['Motif bookkeeping','period ('+Math.floor((r-1)/4)+', '+Math.floor((c-1)/2)+'), local row '+((r-1)%4+1)+', column '+((c-1)%2+1)],
      ['Easy-axis vector e',vec(s.axes[selected])],['Initial spin s₀',vec(s.spins[selected])],
      ['Receiver scale h',s.scale[selected]+' · raw code mask'],['DMI masks 1 / 2 / 3',s.dmiMasks[selected].join(' / ')]];
    $('site-details').innerHTML='<dl>'+entries.map(([k,v])=>'<div><dt>'+esc(k)+'</dt><dd>'+esc(v)+'</dd></div>').join('')+'</dl>';
    const e=s.axes[selected], v=s.spins[selected];
    $('spin-inset').innerHTML='<defs><marker id="inset-arrow" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><path d="M0 0L7 3.5L0 7" fill="#18313b"/></marker></defs><path d="M40 85H195M110 130V12" stroke="#9ba9a6"/><text x="199" y="90" font-size="12">y</text><text x="115" y="16" font-size="12">z</text><circle cx="110" cy="85" r="2" fill="#18313b"/>'+(g==='empty'?'':
      '<line x1="'+(110-45*e[1])+'" y1="'+(85+45*e[2])+'" x2="'+(110+45*e[1])+'" y2="'+(85-45*e[2])+'" stroke="'+colors[g]+'" stroke-width="3" stroke-dasharray="5 4"/><line x1="110" y1="85" x2="'+(110+48*v[1])+'" y2="'+(85-48*v[2])+'" stroke="#18313b" stroke-width="2" marker-end="url(#inset-arrow)"/>')+'<text x="10" y="147" font-size="11" fill="#53696f">solid arrow: initial spin · dashed: easy axis</text>';
  }
  function renderSVG() {
    const svg=$('lattice'), reads=activeReads(), selectedRead=reads.find(x=>x.q===readQ);
    const narrow=window.innerWidth<=720;
    svg.setAttribute('viewBox',narrow?'160 0 540 540':'0 0 860 540');
    const scale=narrow?Math.min(440/(s.L+.5),360/((s.W-1)*Math.sqrt(3)/2+1.8)):Math.min(660/(s.L+1),390/((s.W+1)*Math.sqrt(3)/2));
    const baseX=430-(s.L-.5)*scale/2, baseY=270-(s.W-1)*Math.sqrt(3)/2*scale/2;
    const pos=(r,c)=>{const p=M.position(r,c);return [baseX+p[0]*scale,baseY+p[1]*scale];};
    const point=i=>pos(...M.rc(s,i));
    const [tx,ty]=point(selected), zoom=Number($('zoom').value)/100;
    let out='<defs><marker id="flow" markerWidth="7" markerHeight="7" refX="10" refY="3.5" orient="auto"><path d="M0 0L7 3.5L0 7" fill="context-stroke"/></marker></defs>';
    out+='<text x="'+(narrow?175:24)+'" y="28" font-size="'+(narrow?15:12)+'" fill="#53696f">columns → · rows ↓ · schematic spacing</text>';
    out+='<g transform="translate('+(430+pan[0])+' '+(270+pan[1])+') scale('+zoom+') translate(-430 -270)">';
    if(isCell && $('frames').checked) {
      for(let r=1;r<s.W;r+=4) for(let c=1;c<=s.L;c+=2) {
        const x=baseX+(c-1-.25)*scale,y=baseY+((r-1)*Math.sqrt(3)/2-.25)*scale;
        out+='<rect x="'+x+'" y="'+y+'" width="'+(Math.min(2,s.L-c+1)*scale)+'" height="'+(4*Math.sqrt(3)/2*scale)+'" fill="none" stroke="#9da9a1" stroke-dasharray="5 5"/><text x="'+(x+5)+'" y="'+(y+14)+'" font-size="10" fill="#71827a">'+Math.floor((r-1)/4)+','+Math.floor((c-1)/2)+'</text>';
      }
    }
    // Context connections are the open-boundary union, never a reciprocal-energy claim.
    const seen=new Set();
    for(let i=0;i<s.occupancy.length;i++) if(s.occupancy[i]===1) {
      for(const read of M.reads(s,i,'exchange',1).filter(r=>r.active)) {
        const key=[Math.min(i,read.j),Math.max(i,read.j)].join(':'); if(seen.has(key))continue;seen.add(key);
        const a=point(i),b=point(read.j);
        out+='<line x1="'+a[0]+'" y1="'+a[1]+'" x2="'+b[0]+'" y2="'+b[1]+'" stroke="'+(isCell?'#c0cbc5':'#e0e4db')+'" stroke-width="1.5"/>';
      }
    }
    if(!isCell) for(const read of reads.filter(r=>r.active || r.q===readQ)) {
      const p=pos(read.rr,read.cc), chosen=read.q===readQ;
      const color=chosen?'#b76b10':!read.active?'#939c98':read.weight>0?'#226a9b':'#b04642';
      out+='<g data-read="'+read.q+'" tabindex="0" role="button" aria-label="Inspect '+read.name+' read"><line x1="'+p[0]+'" y1="'+p[1]+'" x2="'+tx+'" y2="'+ty+'" stroke="'+color+'" stroke-width="'+(chosen?4:2.5)+'" '+(read.wrapped||!read.active?'stroke-dasharray="6 4" ':'')+'marker-end="url(#flow)"/><text x="'+((tx+p[0])/2+7)+'" y="'+((ty+p[1])/2-7)+'" font-size="12" fill="'+color+'">'+(read.q+1)+'</text>';
      if(read.outside) out+='<circle cx="'+p[0]+'" cy="'+p[1]+'" r="10" fill="#fffdf8" stroke="'+color+'" stroke-dasharray="3 3"/><text x="'+p[0]+'" y="'+(p[1]-15)+'" text-anchor="middle" font-size="10" fill="'+color+'">'+(read.wrapped?'image → '+pair(read.j):'outside')+'</text>';
      out+='</g>';
    }
    for(let i=0;i<s.occupancy.length;i++) {
      const g=M.group(s,i), empty=g==='empty';
      if(empty && !$('padding').checked && i!==selected)continue;
      const [x,y]=point(i), radius=Math.min(9,scale*.16), neighbor=selectedRead?.j===i;
      out+='<g class="site" data-site="'+i+'" role="button" tabindex="0" aria-label="'+(empty?'Empty slot ':'Mn class '+g+' ')+pair(i)+'"><title>'+pair(i)+' · '+(empty?'empty':'class '+g)+'</title>';
      if(i===selected||neighbor)out+='<circle cx="'+x+'" cy="'+y+'" r="'+(radius+5)+'" fill="none" stroke="'+(i===selected?'#b76b10':'#226a9b')+'" stroke-width="3"/>';
      out+='<circle cx="'+x+'" cy="'+y+'" r="'+radius+'" fill="'+(empty?'#fffdf8':colors[g])+'" stroke="'+colors[g]+'" stroke-width="1.5" '+(empty?'stroke-dasharray="3 2"':'')+'/>';
      if(!empty && $('axes').checked) {
        const e=s.axes[i],d=Math.min(24,scale*.35);
        out+='<line x1="'+(x-d*e[1])+'" y1="'+(y+d*e[2])+'" x2="'+(x+d*e[1])+'" y2="'+(y-d*e[2])+'" stroke="'+colors[g]+'" stroke-width="2" stroke-dasharray="4 3"/>';
      }
      if(!empty && $('spins').checked) {
        const v=s.spins[i],d=Math.min(26,scale*.4);
        out+='<line x1="'+x+'" y1="'+y+'" x2="'+(x+d*v[1])+'" y2="'+(y-d*v[2])+'" stroke="#18313b" stroke-width="2" marker-end="url(#flow)"/>';
      }
      if($('labels').checked)out+='<text x="'+x+'" y="'+(y+radius+16)+'" text-anchor="middle" font-size="'+(narrow?14:s.W>9?9:11)+'" fill="#324d55">'+M.rc(s,i).join(',')+'</text>';
      out+='</g>';
    }
    out+='</g>';svg.innerHTML=out;
  }
  function renderArray() {
    const reads=activeReads(), selectedRead=reads.find(r=>r.q===readQ);
    const neighbors=new Set(selectedRead ? [selectedRead.j] : isCell?[]:reads.filter(r=>r.active).map(r=>r.j));
    let html='<table class="array" aria-label="MATLAB occupation array"><thead><tr><th>r / c</th>';
    for(let c=1;c<=s.L;c++)html+='<th scope="col">'+c+'</th>';
    html+='</tr></thead><tbody>';
    for(let r=1;r<=s.W;r++) {
      html+='<tr><th scope="row">'+r+'</th>';
      for(let c=1;c<=s.L;c++) {
        const i=M.index(s,r,c),g=M.group(s,i),empty=g==='empty';
        const classes=[i===selected?'selected':'',neighbors.has(i)?'neighbor':'',empty?'empty':'',empty&&!$('padding').checked&&i!==selected?'hidden-padding':''].join(' ');
        html+='<td><button type="button" data-site="'+i+'" class="'+classes+'" style="'+(empty?'':'color:'+colors[g]+';background:'+colors[g]+'14')+'" aria-label="'+pair(i)+', '+(empty?'empty':'class '+g)+'" aria-pressed="'+(i===selected)+'">'+(empty?'·':g)+'</button></td>';
      }html+='</tr>';
    }$('array-view').innerHTML=html+'</tbody></table>';
  }
  function renderTable() {
    const t=term(), reads=activeReads(), p=data.parameters;
    $('reads-table').querySelector('thead').innerHTML=t==='anisotropy'?'<tr><th>Entry</th><th>e = (x,y,z)</th><th>K / μ (T)</th><th>Status</th></tr>':'<tr><th>Read</th><th>Source offset<br>Δrow, Δcol</th><th>Source j</th><th>Raw mask</th><th>h × mask</th><th>Energy-scale coefficient (J)</th><th>Field coefficient (T)</th><th>Reverse weight</th><th>Status</th></tr>';
    if(t==='anisotropy') {
      $('reads-table').querySelector('tbody').innerHTML='<tr><td>'+pair(selected)+'</td><td>'+vec(s.axes[selected])+'</td><td>'+sci(p.K/p.mu)+'</td><td>'+(s.occupancy[selected]===1?'onsite':'empty receiver; excluded')+'</td></tr>';
      $('table-note').textContent='Anisotropy is onsite: no neighbor, shift, or bond multiplicity. Empty slots have zero axes and no physical anisotropy response.';
      $('read-detail').innerHTML='<p class="small">Onsite term · select another array entry to inspect its easy axis.</p>';return;
    }
    $('table-note').textContent='All eight stencil directions are listed, including cuts and empty slots. Coefficients are raw algebraic values before exclusion; only “active” rows contribute to this occupied receiver. Reverse weight is the summed active coefficient for i read by j, in units of '+(t==='exchange'?'Jintra':'Dsim')+'.';
    $('reads-table').querySelector('tbody').innerHTML=reads.map(read=>{
      const energy=read.weight*(t==='exchange'?p.J:p.D), field=energy/p.mu*(t==='exchange'?-1:1);
      return '<tr class="'+(readQ===read.q?'selected':'')+'"><td><button type="button" data-read="'+read.q+'" aria-pressed="'+(readQ===read.q)+'">'+(read.q+1)+' / '+read.name+'</button></td><td>('+signed(read.dr)+', '+signed(read.dc)+')</td><td>'+pair(read.j)+(read.wrapped?'<br><span class="tag">wrapped</span>':'')+'</td><td>'+signed(read.raw)+'</td><td>'+signed(read.weight)+'</td><td>'+sci(energy)+'</td><td>'+sci(field)+(t==='dmi'?'<br><span class="small">By ← sz; Bz ← sy has opposite sign</span>':'')+'</td><td>'+signed(M.reverse(s,selected,read,t,bc()))+'</td><td class="read-status">'+read.status+'</td></tr>';
    }).join('');
    const read=reads.find(x=>x.q===readQ);
    $('read-detail').innerHTML=read?'<h3>Read '+(read.q+1)+' / '+read.name+'</h3><p class="small">Requested index ('+read.rr+', '+read.cc+') → '+pair(read.j)+'.<br>Array-wrap translation: ('+read.wrap.join(', ')+') periods.<br>'+read.status+'.'+(read.wrapped?' The ghost marker shows the requested image; the array highlights its wrapped endpoint.':'')+'</p>':'<p class="small">Select a neighbor row or a numbered link to inspect its path.</p>';
  }
  function renderFormula() {
    const t=term(), probe=s.probes.find(p=>p.bc===bc());
    $('boundary-note').textContent=bc()===1?'Open mode: off-array reads are zeroed. The occupation mask also removes selected last-column sites.':'Legacy bc = 0: shifts wrap around the array, but occupation and exchange masks retain edge cuts. This is not a verified periodic kagome cell; exchange and DMI seam behavior differ.';
    $('diagram-heading').textContent=t==='anisotropy'?'01 / Onsite axis':'01 / Selected '+(t==='dmi'?'signed DMI':'exchange')+' reads';
    if(t==='exchange') {
      $('formula').innerHTML='<div class="formula">B<sub>i</sub><sup>ex</sup> = −(h<sub>i</sub> / μ<sub>i</sub>) Σ<sub>q</sub> A<sub>i,q</sub> s<sub>j(i,q)</sub></div><p>Each coefficient belongs to the receiving site. Positive raw weights enter with a minus sign in the field. The table lists both directions separately.</p>';
      $('term-note').textContent='Jintra = '+sci(data.parameters.J)+' J; μ = '+sci(data.parameters.mu)+' J/T. This patch has '+probe.unequalPairs+' unequal reverse pairs (max difference '+probe.maxExchangeAsymmetry+' × Jintra). A reciprocal pair-energy interpretation remains unresolved; no extra half factor is introduced.';
    }else if(t==='dmi') {
      $('formula').innerHTML='<div class="formula">B<sub>x</sub> = 0<br>B<sub>y,i</sub> = (h<sub>i</sub>D / μ<sub>i</sub>) Σ<sub>q</sub> c<sub>i,q</sub>s<sub>z,j</sub><br>B<sub>z,i</sub> = −(h<sub>i</sub>D / μ<sub>i</sub>) Σ<sub>q</sub> c<sub>i,q</sub>s<sub>y,j</sub></div><p>Blue and red distinguish positive and negative c coefficients. These are signed component reads; they do not establish physical DMI-vector directions.</p>';
      $('term-note').textContent='D = '+sci(data.parameters.D)+' J. The maximum |C + Cᵀ| of the occupied-site normalized operator is '+probe.dmiAntisymmetryError+'. This diagnostic includes receiver scaling. Physical bond orientation, energy sign and edge normalization remain unresolved.';
    }else {
      $('formula').innerHTML='<div class="formula">B<sub>i</sub><sup>an</sup> = (K / μ<sub>i</sub>)(e<sub>i</sub> · s<sub>i</sub>) e<sub>i</sub></div><p>The axis is a unit vector at occupied entries. This current field follows algebraically from E<sub>i</sub> = −½K(e<sub>i</sub> · s<sub>i</sub>)² under B = −(1/μ)∂E/∂s.</p>';
      $('term-note').textContent='K = '+sci(data.parameters.K)+' J. The factor relative to the published single-site anisotropy remains unresolved; this page preserves the implemented coefficient.';
    }
  }
  function excerpt(file,start,end,highlight) {
    const lines=data.sources[file].split(/\r?\n/),a=lines.findIndex(l=>l.includes(start));
    let b=lines.findIndex((l,i)=>i>a&&l.includes(end)); if(b<0)b=lines.length;
    if(a<0)return '';
    const body=lines.slice(a,b).map((l,i)=>{
      const value=esc(String(a+i+1).padStart(3)+'  '+l);
      return highlight && l.includes(highlight)?'<mark>'+value+'</mark>':value;
    }).join('\n');
    return '<h3><a href="../'+file+'.m">'+file+'.m</a> · lines '+(a+1)+'–'+b+'</h3><pre>'+body+'</pre>';
  }
  function sourceExcerpts() {
    const q=readQ>=0?data.shiftNames[readQ]:null;
    if(isCell) {
      $('source-excerpts').innerHTML=excerpt('distrib','%% Mn3Sn','%% hex')+excerpt('systemgeneration','            phi_=90/180*pi;','clear ctL ctW');
    }else if(term()==='exchange') {
      $('source-excerpts').innerHTML=excerpt('integrate_llg','AsimnextL(2:2','for ctL=1:natomL',q)+excerpt('distrib','% mid/2','%% hdmi')+excerpt('integrate_llg','        mmxnextL=circshift','        % hk_atom',q);
    }else if(term()==='dmi') {
      $('source-excerpts').innerHTML=excerpt('distrib','%% hdmi','%% hani')+excerpt('integrate_llg','        mmxnextL=circshift','        hex_x=',q)+excerpt('integrate_llg','        hdmi_x=','        if thermalenable',q);
    }else {
      $('source-excerpts').innerHTML=excerpt('distrib','%% hani','if load_fixed_atom_distrib')+excerpt('integrate_llg','        hani_x=','        % Renjie');
    }
  }
  function provenance() {
    $('provenance').innerHTML='<p>Generated by <a href="export_visualization_data.m">export_visualization_data.m</a> using MATLAB '+esc(data.matlab)+'. Every occupied Cartesian impulse was checked for exchange, DMI, and anisotropy in both boundary modes at normalized tolerance '+data.tolerance+'. This validates source correspondence, not physical correctness.</p><p>Refresh and check instructions: <a href="../README.md#interactive-structure-and-interactions">README</a>. SHA-256 source hashes:</p><ul>'+Object.entries(data.hashes).map(([k,v])=>'<li>'+k+'.m: <code>'+v+'</code></li>').join('')+'</ul>';
  }
  $('preset').addEventListener('change',()=>setSize(...$('preset').value.split('x').map(Number)));
  if(isCell) for(const id of ['repeat-w','repeat-l']) $(id).addEventListener('change',()=>{
    const W=Number($('repeat-w').value),L=Number($('repeat-l').value);setSize(W,W!==5&&L===5?6:L);
  });
  $('target').addEventListener('change',()=>selectSite(Number($('target').value)));
  if(!isCell) for(const id of ['term','boundary']) $(id).addEventListener('change',()=>{readQ=-1;render();});
  for(const id of ['labels','padding','frames','spins','axes']) if($(id))$(id).addEventListener('change',()=>{renderSVG();renderArray();});
  $('zoom').addEventListener('input',renderSVG);
  $('reset').addEventListener('click',()=>{pan=[0,0];$('zoom').value=100;renderSVG();});
  document.addEventListener('click',event=>{
    if(moved){moved=false;return;}
    const site=event.target.closest('[data-site]'),read=event.target.closest('[data-read]'),bond=event.target.closest('[data-exchange-pair]');
    if(bond){
      const map=bond.closest('svg'),value=bond.dataset.exchangePair;
      selectExchangePair(value);
      (map?$('unequal-map'):$('unequal-table')).querySelector('[data-exchange-pair="'+value+'"]').focus({preventScroll:true});
    }
    else if(site)selectSite(Number(site.dataset.site));
    else if(read){readQ=Number(read.dataset.read);render();}
  });
  if(!isCell)$('unequal-map').addEventListener('keydown',event=>{
    const bond=event.target.closest('[data-exchange-pair]');
    if(bond&&(event.key==='Enter'||event.key===' ')){event.preventDefault();bond.dispatchEvent(new MouseEvent('click',{bubbles:true}));}
  });
  $('lattice').addEventListener('keydown',event=>{
    if(event.key==='Enter'||event.key===' '){const item=event.target.closest('[data-site],[data-read]');if(item){event.preventDefault();item.dispatchEvent(new MouseEvent('click',{bubbles:true}));return;}}
    if(['+','=','-','0','ArrowLeft','ArrowRight','ArrowUp','ArrowDown'].includes(event.key)){
      event.preventDefault();
      if(event.key==='0'){pan=[0,0];$('zoom').value=100;}
      else if(event.key.startsWith('Arrow')){pan[0]+=event.key==='ArrowLeft'?-20:event.key==='ArrowRight'?20:0;pan[1]+=event.key==='ArrowUp'?-20:event.key==='ArrowDown'?20:0;}
      else $('zoom').value=Math.max(70,Math.min(220,Number($('zoom').value)+(event.key==='-'?-10:10)));
      renderSVG();
    }
  });
  $('lattice').addEventListener('pointerdown',event=>{drag=[event.clientX,event.clientY,...pan];moved=false;});
  window.addEventListener('pointermove',event=>{
    if(!drag)return;const dx=event.clientX-drag[0],dy=event.clientY-drag[1];
    if(Math.abs(dx)+Math.abs(dy)>4){moved=true;const rect=$('lattice').getBoundingClientRect();const ratio=Math.min(rect.width/(window.innerWidth<=720?540:860),rect.height/540);pan=[drag[2]+dx/ratio,drag[3]+dy/ratio];renderSVG();}
  });
  window.addEventListener('pointerup',()=>{drag=null;setTimeout(()=>{moved=false;},0);});
  window.addEventListener('pointercancel',()=>{drag=null;moved=false;});
  window.addEventListener('resize',renderSVG);
  syncControls();provenance();render();
})();
