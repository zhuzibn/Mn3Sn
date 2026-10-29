/* Run from the repository root: node docs/check_visualization_data.cjs */
'use strict';
const fs=require('node:fs'), path=require('node:path'), vm=require('node:vm');
const assert=require('node:assert/strict'), crypto=require('node:crypto');
const root=path.resolve(__dirname,'..'), context={window:{}};
for(const file of ['mn3sn-data.js','mn3sn-model.js'])vm.runInNewContext(fs.readFileSync(path.join(__dirname,file),'utf8'),context);
const M=context.window.Mn3Sn, data=M.data;
for(const [key,hash] of Object.entries(data.hashes)) {
  const bytes=fs.readFileSync(path.join(root,key+'.m'));
  assert.equal(crypto.createHash('sha256').update(bytes).digest('hex'),hash,'Stale snapshot: '+key+'.m');
  assert.equal(bytes.toString('utf8'),data.sources[key],'Stale source excerpt: '+key+'.m');
}
let comparisons=0,activeWraps=0,emptySources=0,emptyReceivers=0,cuts=0,openEdges=0;
for(const s of data.snapshots) {
  const n=s.W*s.L;
  for(let i=0;i<n;i++) {
    const [r,c]=M.rc(s,i); assert.equal(M.index(s,r,c),i);
    const empty=(r%4===1&&c%2===0)||(r%4===3&&c%2===1)||(r%2===1&&c===s.L);
    assert.equal(s.occupancy[i],empty?3:1,'Occupation at '+[s.W,s.L,r,c]);
    if(empty) {assert(s.spins[i].every(x=>x===0));assert(s.axes[i].every(x=>x===0));}
    else {
      const group=M.group(s,i), angle={A:30,B:150,C:270}[group]*Math.PI/180;
      assert(Math.abs(s.axes[i][1]-Math.sin(angle))<1e-12);
      assert(Math.abs(s.axes[i][2]-Math.cos(angle))<1e-12);
      assert(Math.abs(s.spins[i].reduce((sum,x)=>sum+x*x,0)-1)<1e-12);
    }
  }
  for(const probe of s.probes) {
    const expected=new Map(probe.operator.map(row=>[row[0]+':'+row[1],row.slice(2)]));
    const actual=[new Float64Array(n*n),new Float64Array(n*n)];
    for(const [t,term] of ['exchange','dmi'].entries())for(let i=0;i<n;i++) {
      const reads=M.reads(s,i,term,probe.bc); assert.equal(reads.length,8);
      for(const read of reads) {
        if(read.status==='coefficient cut')cuts++;
        if(read.status==='open edge')openEdges++;
        if(read.status==='empty receiver')emptyReceivers++;
        if(read.status==='empty source')emptySources++;
        if(read.active){actual[t][i*n+read.j]+=read.weight;if(read.wrapped)activeWraps++;}
        if(probe.bc===1&&read.outside)assert.equal(read.active,false);
        if(read.active)assert.equal(s.occupancy[read.j],1);
      }
    }
    for(let i=0;i<n;i++)for(let j=0;j<n;j++)for(let t=0;t<2;t++) {
      const wanted=expected.get(i+':'+j)?.[t]||0;
      assert(Math.abs(actual[t][i*n+j]-wanted)<1e-12,'Impulse mismatch '+[s.W,s.L,probe.bc,i,j,t]);comparisons++;
      if(wanted===0)assert.equal(actual[t][i*n+j],0);
    }
    let pairs=0,mismatches=0;
    for(let i=0;i<n;i++)for(let j=i+1;j<n;j++) {
      if(actual[0][i*n+j]||actual[0][j*n+i])pairs++;
      if(Math.abs(actual[0][i*n+j]-actual[0][j*n+i])>1e-12)mismatches++;
    }
    assert.equal(pairs,probe.exchangePairs);assert.equal(mismatches,probe.unequalPairs);
    const overview=M.exchangePairs(s,probe.bc);
    assert.equal(overview.length,pairs,'Overview pair count');
    assert.equal(overview.filter(p=>Math.abs(p.intoA-p.intoB)>data.tolerance).length,mismatches,'Overview unequal count');
    assert.equal(new Set(overview.map(p=>p.a+':'+p.b)).size,overview.length,'Duplicate overview pairs');
    for(const p of overview) {
      assert(p.a<p.b,'Overview endpoint order');
      assert.equal(s.occupancy[p.a],1);assert.equal(s.occupancy[p.b],1);
      assert(Math.abs(p.intoA-(expected.get(p.a+':'+p.b)?.[0]||0))<1e-12,'Overview weight into first endpoint');
      assert(Math.abs(p.intoB-(expected.get(p.b+':'+p.a)?.[0]||0))<1e-12,'Overview weight into second endpoint');
    }
    for(let i=0;i<n;i++)for(const term of ['exchange','dmi'])for(const read of M.reads(s,i,term,probe.bc)) {
      const reverse=read.j===null?0:actual[term==='exchange'?0:1][read.j*n+i];
      assert.equal(M.reverse(s,i,read,term,probe.bc),reverse);
    }
  }
}
assert(activeWraps>0&&emptySources>0&&emptyReceivers>0&&cuts>0&&openEdges>0,'Missing boundary/status coverage');
const notebook=M.snapshot(5,6),open=notebook.probes.find(p=>p.bc===1);
assert.equal(notebook.occupancy.filter(x=>x===1).length,20);
assert.equal(open.exchangePairs,30);assert.equal(open.unequalPairs,12);
console.log('PASS: source hashes/excerpts, '+data.snapshots.length+' sizes × 2 boundaries, '+comparisons+' operator entries against independent MATLAB impulses.');
console.log({activeWraps,emptySources,emptyReceivers,cuts,openEdges});
