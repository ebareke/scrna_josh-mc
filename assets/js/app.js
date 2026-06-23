/* ============================================================================
 * E9.5 Eftud2 / Sf3b4 / Control scRNA-seq dashboard
 * Author: Eric Bareke
 * ==========================================================================*/
'use strict';
const D = window.DASH_DATA;
const COND = ['Control', 'Eftud2', 'Sf3b4'];
const COL  = {Control:'#3b4992', Eftud2:'#ee0000', Sf3b4:'#008b45'};
const COLM = {Control:'#6a78c9', Eftud2:'#ff5a4d', Sf3b4:'#2fb877'}; // lighter
const PALETTE = ['#5b8def','#33c7b5','#e0533d','#f2b134','#9b6cf0','#3fd49b','#ef6fb0',
                 '#6ad1e0','#c0d24a','#f08a5d','#8a93c2','#d35d6e','#5cc99a','#b388eb'];
let DS = 'ALL_CELLS';
const $ = (s,r=document)=>r.querySelector(s);
const el = (t,c,h)=>{const e=document.createElement(t);if(c)e.className=c;if(h!=null)e.innerHTML=h;return e;};
const pretty = s => (s||'').replace(/_/g,' ');
const cmpLabel = c => c.replace('_vs_',' vs ').replace(/_/g,' ');
const fmtP = p => p==null?'—':(p<1e-4?p.toExponential(1):p.toFixed(p<0.01?4:3));
const fmtN = (n,d=2)=>n==null?'—':(+n).toFixed(d);

/* ---------------- navigation ---------------- */
const SECTIONS = [
  ['overview','Overview', icon('grid')],
  ['composition','Cell-type composition', icon('chart')],
  ['de','Differential expression', icon('volcano')],
  ['pseudobulk','Pseudo-bulk DE', icon('layers')],
  ['enrichment','Functional enrichment', icon('target')],
  ['splicing','Splicing ratios', icon('dna')],
  ['velocity','RNA velocity', icon('flow')],
  ['methods','Methods &amp; data', icon('book')],
];
function icon(n){
  const p={
    grid:'<rect x="3" y="3" width="7" height="7" rx="1"/><rect x="14" y="3" width="7" height="7" rx="1"/><rect x="3" y="14" width="7" height="7" rx="1"/><rect x="14" y="14" width="7" height="7" rx="1"/>',
    chart:'<path d="M4 20V10M10 20V4M16 20v-7M22 20H2"/>',
    volcano:'<path d="M12 3v6M5 21l4-7M19 21l-4-7"/><circle cx="12" cy="11" r="1.5"/><circle cx="8" cy="16" r="1.2"/><circle cx="16" cy="15" r="1.2"/>',
    layers:'<path d="M12 3l9 5-9 5-9-5 9-5zM3 14l9 5 9-5"/>',
    target:'<circle cx="12" cy="12" r="8"/><circle cx="12" cy="12" r="3.5"/>',
    dna:'<path d="M7 3c0 5 10 5 10 10S7 18 7 21M17 3c0 5-10 5-10 10s10 3 10 8"/>',
    flow:'<path d="M4 7h10a3 3 0 010 6H8a3 3 0 000 6h12"/><path d="M16 4l4 3-4 3"/>',
    book:'<path d="M4 5a2 2 0 012-2h12v18H6a2 2 0 01-2-2zM8 7h8M8 11h8"/>',
  }[n]||'';
  return `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">${p}</svg>`;
}
function buildNav(){
  const nav=$('#nav'); nav.innerHTML='';
  nav.appendChild(el('div','grp','Explore'));
  SECTIONS.forEach(([id,label,ic],i)=>{
    const a=el('a',i===0?'active':'',`${ic}<span>${label}</span>`);
    a.href='#'+id; a.dataset.sec=id;
    a.onclick=e=>{e.preventDefault();showSection(id);};
    nav.appendChild(a);
    if(id==='overview'||id==='enrichment') nav.appendChild(el('div','grp',id==='overview'?'Differences':'More'));
  });
}
function showSection(id){
  document.querySelectorAll('nav.side a').forEach(a=>a.classList.toggle('active',a.dataset.sec===id));
  document.querySelectorAll('.section').forEach(s=>s.classList.toggle('active',s.dataset.sec===id));
  window.scrollTo({top:0,behavior:'smooth'});
  setTimeout(()=>window.dispatchEvent(new Event('resize')),60);
}

/* ---------------- lightbox ---------------- */
function lightbox(src,cap){
  $('#lbimg').src=src; $('#lbcap').textContent=cap||''; $('#lb').classList.add('open');
}
$('#lb').onclick=e=>{if(e.target.id!=='lbimg')$('#lb').classList.remove('open');};
$('#lbx').onclick=()=>$('#lb').classList.remove('open');
document.addEventListener('keydown',e=>{if(e.key==='Escape')$('#lb').classList.remove('open');});

function thumb(src,cap){
  if(!src) return null;
  const t=el('div','thumb');
  t.innerHTML=`<img loading="lazy" src="${src}" alt="${cap}"><div class="cap">${cap}</div>`;
  t.onclick=()=>lightbox(src,cap);
  return t;
}
function gallery(items){ // items: [[figKey,caption],...]
  const g=el('div','gallery'); let any=false;
  items.forEach(([k,c])=>{const f=fig(k);if(f){const t=thumb(f,c);if(t){g.appendChild(t);any=true;}}});
  return any?g:el('div','empty','No figures available for this view.');
}
const ds=()=>D.datasets[DS];
const fig=k=>ds().figures&&ds().figures[k];

/* common Plotly layout */
function lay(extra={}){
  return Object.assign({
    paper_bgcolor:'rgba(0,0,0,0)', plot_bgcolor:'rgba(0,0,0,0)',
    font:{color:'#c7d1e6',family:'Inter, sans-serif',size:12},
    margin:{l:60,r:18,t:30,b:60},
    xaxis:{gridcolor:'#243049',zerolinecolor:'#2a3450',linecolor:'#2a3450'},
    yaxis:{gridcolor:'#243049',zerolinecolor:'#2a3450',linecolor:'#2a3450'},
    legend:{bgcolor:'rgba(0,0,0,0)',font:{size:11}},
    hoverlabel:{bgcolor:'#1d2640',bordercolor:'#2a3450',font:{color:'#e8edf6'}},
  },extra);
}
const CFG={displayModeBar:true,responsive:true,displaylogo:false,
  modeBarButtonsToRemove:['lasso2d','select2d','autoScale2d','toggleSpikelines']};

/* ============================================================================
 * SECTION RENDERERS
 * ==========================================================================*/
function card(title,hint){
  const c=el('div','card');
  if(title) c.appendChild(el('h3',null,title));
  if(hint) c.appendChild(el('div','hint',hint));
  return c;
}
function plotInto(card,id,h=360){
  const d=el('div','plot'); d.id=id; d.style.minHeight=h+'px'; card.appendChild(d); return d;
}

/* ---------- OVERVIEW ---------- */
function renderOverview(sec){
  sec.innerHTML='';
  const d=ds();
  const dsName = DS==='ALL_CELLS'?'whole-embryo atlas (all cell types)':'neural-crest compartment (NCC sub-types)';
  sec.appendChild(el('h2',null,`Overview — ${DS==='ALL_CELLS'?'All cell types':'Neural crest cells'}`));
  sec.appendChild(el('p','lead',`A bird's-eye comparison of the ${dsName} across the three E9.5 conditions. Mutants <i>Eftud2</i> and <i>Sf3b4</i> are contrasted against the wild-type <b>Control</b> baseline. Use the left-hand menu to drill into composition, differential expression, enrichment and splicing.`));

  // KPI tiles
  const nCT=d.composition.cell_types.length;
  const totSig=(()=>{let u=0,dn=0;Object.values(d.de||{}).forEach(byct=>Object.values(byct).forEach(o=>{u+=o.up;dn+=o.dn;}));return u+dn;})();
  const spl=d.splicing_global;
  const kpis=el('div','kpis');
  const tiles=[
    ['',nCT,'cell '+(DS==='ALL_CELLS'?'types':'sub-types')+' annotated'],
    ['',Object.keys(d.de||{}).length,'pairwise DE contrasts'],
    ['',totSig.toLocaleString(),'significant DE calls (FDR<0.05)'],
    ['c',fmtN(spl.Control,2),'Control splicing ratio'],
    ['e',fmtN(spl.Eftud2,2),'Eftud2 splicing ratio'],
    ['s',fmtN(spl.Sf3b4,2),'Sf3b4 splicing ratio'],
  ];
  tiles.forEach(([cls,v,l])=>{const k=el('div','kpi');k.innerHTML=`<div class="v ${cls}">${v}</div><div class="l">${l}</div>`;kpis.appendChild(k);});
  sec.appendChild(kpis);

  sec.appendChild(el('div','callout',`<b>Splicing signal.</b> Splicing ratio = unspliced ⁄ (spliced+1) per cell, averaged per condition. Both spliceosome mutants show globally elevated ratios versus Control — consistent with disrupted splicing kinetics — quantified per cell type in the Splicing section.`));

  // condition UMAP triptych
  const c1=card('Annotated UMAP per condition','Cells coloured by annotated cell type. Click any panel to enlarge. Each condition pools two replicates.');
  const tri=el('div','cond-tri');
  COND.forEach(cn=>{const t=thumb(fig('umap_'+cn),cn+' — annotated UMAP');if(t)tri.appendChild(t);});
  c1.appendChild(tri.children.length?tri:el('div','empty','UMAPs unavailable.'));
  sec.appendChild(c1);

  // composition + splicing quick charts
  const g=el('div','grid2');
  const cc=card('Cell-type composition shift','Stacked % of cells per condition.');
  plotInto(cc,'ov-comp',380); g.appendChild(cc);
  const sc=card('Global splicing ratio by condition','n-weighted mean of per-cell unspliced/spliced ratio.');
  plotInto(sc,'ov-spl',380); g.appendChild(sc);
  sec.appendChild(g);

  drawCompositionStacked('ov-comp');
  drawSplicingGlobal('ov-spl');
}

/* ---------- COMPOSITION ---------- */
function renderComposition(sec){
  sec.innerHTML='';
  const d=ds();
  sec.appendChild(el('h2',null,'Cell-type composition'));
  sec.appendChild(el('p','lead','How the cellular make-up of the embryo shifts between Control and the two spliceosome mutants. Compositional differences are tested with a CLR linear-mixed model (replicate-aware omnibus + pairwise) and per-cell-type Fisher exact tests (odds ratios).'));

  // view switch
  const ctrl=el('div','controls');
  ctrl.innerHTML=`<label>View</label>`;
  const seg=el('div','seg');
  [['stack','Stacked %'],['group','Grouped bars'],['delta','Mutant − Control']].forEach(([v,l],i)=>{
    const b=el('button',i===0?'active':'',l);b.dataset.v=v;
    b.onclick=()=>{seg.querySelectorAll('button').forEach(x=>x.classList.remove('active'));b.classList.add('active');drawCompositionMain('comp-main',v);};
    seg.appendChild(b);
  });
  ctrl.appendChild(seg);
  const mc=card('Composition across conditions');
  mc.appendChild(ctrl);
  plotInto(mc,'comp-main',440);
  sec.appendChild(mc);

  // stats: fisher OR + CLR pairwise
  const g=el('div','grid2');
  const fc=card('Per-cell-type enrichment — Fisher odds ratios','log₂ OR per cell type for each mutant contrast (>0 enriched in first group). Filled = FDR<0.05.');
  plotInto(fc,'comp-fisher',430); g.appendChild(fc);
  const clr=card('CLR mixed-model effects (replicate-aware)','Δ CLR abundance per pairwise contrast; sized by significance.');
  plotInto(clr,'comp-clr',430); g.appendChild(clr);
  sec.appendChild(g);

  // omnibus table + chisq
  const oc=card('Compositional significance (CLR-LME omnibus)',`Likelihood-ratio test per cell type across conditions. ${d.chisq&&d.chisq[0]?'Pooled χ² test: χ²='+(+d.chisq[0].X2).toFixed(0)+', df='+d.chisq[0].df+', p'+(+d.chisq[0].pvalue===0?'<1e-300':'='+d.chisq[0].pvalue)+'.':''}`);
  oc.appendChild(omnibusTable(d.clr_omnibus));
  sec.appendChild(oc);

  // existing figures
  const fc2=card('Reference composition figures','Publication-style plots generated by the analysis pipeline.');
  fc2.appendChild(gallery([
    ['conditions_composition_ranked_bar','Ranked composition by condition'],
    ['conditions_composition_dotplot','Composition dot-plot'],
    ['04_lollipop_both_mutants','Mutant vs Control lollipop'],
    ['composition_OR_heatmap','Odds-ratio heatmap'],
    ['composition_clr_lme_forest','CLR-LME forest plot'],
    ['05_bubble_composition_signgene','Composition × significant-gene bubble'],
  ]));
  sec.appendChild(fc2);

  drawCompositionMain('comp-main','stack');
  drawFisher('comp-fisher');
  drawCLR('comp-clr');
}

function drawCompositionStacked(id){
  const c=ds().composition, cts=c.cell_types;
  const traces=cts.map((ct,i)=>({
    x:COND, y:COND.map(cn=>c[cn][i]), name:pretty(ct), type:'bar',
    marker:{color:PALETTE[i%PALETTE.length]}, hovertemplate:`%{x}<br>${pretty(ct)}: %{y:.1f}%<extra></extra>`,
  }));
  Plotly.newPlot(id,traces,lay({barmode:'stack',yaxis:{title:'% of cells',gridcolor:'#243049'},
    legend:{font:{size:10},orientation:'v'},margin:{l:55,r:10,t:10,b:40}}),CFG);
}
function drawSplicingGlobal(id){
  const g=ds().splicing_global;
  Plotly.newPlot(id,[{x:COND,y:COND.map(c=>g[c]),type:'bar',
    marker:{color:COND.map(c=>COL[c])},text:COND.map(c=>g[c]==null?'':g[c].toFixed(2)),
    textposition:'outside',hovertemplate:'%{x}: %{y:.3f}<extra></extra>'}],
    lay({yaxis:{title:'unspliced / (spliced+1)',gridcolor:'#243049'},margin:{l:60,r:10,t:20,b:40}}),CFG);
}
function drawCompositionMain(id,view){
  const c=ds().composition, cts=c.cell_types;
  let traces,layout;
  if(view==='stack'){
    traces=cts.map((ct,i)=>({x:COND,y:COND.map(cn=>c[cn][i]),name:pretty(ct),type:'bar',
      marker:{color:PALETTE[i%PALETTE.length]},hovertemplate:`%{x} · ${pretty(ct)}: %{y:.1f}%<extra></extra>`}));
    layout=lay({barmode:'stack',yaxis:{title:'% of cells'},margin:{l:60,r:10,t:10,b:50}});
  }else if(view==='group'){
    traces=COND.map(cn=>({x:cts.map(pretty),y:cts.map((ct,i)=>c[cn][i]),name:cn,type:'bar',
      marker:{color:COL[cn]},hovertemplate:`${cn} · %{x}: %{y:.1f}%<extra></extra>`}));
    layout=lay({barmode:'group',yaxis:{title:'% of cells'},xaxis:{tickangle:-35},margin:{l:60,r:10,t:10,b:120}});
  }else{ // delta vs control
    traces=['Eftud2','Sf3b4'].map(cn=>({x:cts.map(pretty),
      y:cts.map((ct,i)=>c[cn][i]-c.Control[i]),name:cn+' − Control',type:'bar',
      marker:{color:COL[cn]},hovertemplate:`${cn}−Control · %{x}: %{y:+.1f} pp<extra></extra>`}));
    layout=lay({barmode:'group',yaxis:{title:'Δ percentage points'},xaxis:{tickangle:-35},margin:{l:60,r:10,t:10,b:120}});
  }
  Plotly.newPlot(id,traces,layout,CFG);
}
function drawFisher(id){
  const f=ds().fisher; if(!f||!f.length){return Plotly.purge(id);}
  const comps=[...new Set(f.map(r=>r.comparison))];
  const traces=comps.map((cmp,ci)=>{
    const rows=f.filter(r=>r.comparison===cmp&&r.log2OR!=null);
    return {x:rows.map(r=>r.log2OR),y:rows.map(r=>pretty(r.cell_type)),
      name:cmpLabel(cmp),mode:'markers',type:'scatter',
      marker:{size:rows.map(r=>r.padj!=null&&r.padj<0.05?13:8),
        color:PALETTE[ci],line:{width:rows.map(r=>r.padj!=null&&r.padj<0.05?2:0),color:'#fff'},
        opacity:rows.map(r=>r.padj!=null&&r.padj<0.05?1:.45)},
      text:rows.map(r=>`${cmpLabel(cmp)}<br>${pretty(r.cell_type)}<br>log₂OR=${fmtN(r.log2OR)} · FDR=${fmtP(r.padj)}`),
      hovertemplate:'%{text}<extra></extra>'};
  });
  Plotly.newPlot(id,traces,lay({xaxis:{title:'log₂ odds ratio',zeroline:true,zerolinecolor:'#52608a'},
    margin:{l:160,r:10,t:10,b:45},legend:{orientation:'h',y:1.08}}),CFG);
}
function drawCLR(id){
  const c=ds().clr_pairwise; if(!c||!c.length){return Plotly.purge(id);}
  const comps=[...new Set(c.map(r=>r.comparison))];
  const traces=comps.map((cmp,ci)=>{
    const rows=c.filter(r=>r.comparison===cmp&&r.delta_clr!=null);
    return {x:rows.map(r=>r.delta_clr),y:rows.map(r=>pretty(r.cell_type)),name:cmpLabel(cmp),
      mode:'markers',type:'scatter',
      marker:{size:rows.map(r=>r.padj!=null&&r.padj<0.05?14:8),color:PALETTE[ci],
        opacity:rows.map(r=>r.padj!=null&&r.padj<0.05?1:.45),line:{width:rows.map(r=>r.padj!=null&&r.padj<0.05?2:0),color:'#fff'}},
      text:rows.map(r=>`${cmpLabel(cmp)}<br>${pretty(r.cell_type)}<br>ΔCLR=${fmtN(r.delta_clr)} · FDR=${fmtP(r.padj)}`),
      hovertemplate:'%{text}<extra></extra>'};
  });
  Plotly.newPlot(id,traces,lay({xaxis:{title:'Δ CLR abundance',zeroline:true,zerolinecolor:'#52608a'},
    margin:{l:160,r:10,t:10,b:45},legend:{orientation:'h',y:1.08}}),CFG);
}
function omnibusTable(rows){
  const s=el('div','tbl-scroll');
  const t=el('table','tbl');
  t.innerHTML='<thead><tr><th>Cell type</th><th class="num">LRT χ²</th><th class="num">df</th><th class="num">p-value</th><th class="num">FDR</th><th>Sig.</th></tr></thead>';
  const tb=el('tbody');
  (rows||[]).slice().sort((a,b)=>(a.padj??9)-(b.padj??9)).forEach(r=>{
    const sig=r.padj!=null&&r.padj<0.05;
    const tr=el('tr',null,`<td class="gene-name">${pretty(r.cell_type)}</td>
      <td class="num">${fmtN(r.chi2,1)}</td><td class="num">${r.padj!=null?'':''}${rowDf(r)}</td>
      <td class="num">${fmtP(r.pvalue)}</td><td class="num">${fmtP(r.padj)}</td>
      <td>${sig?'<span class="tag up">FDR&lt;0.05</span>':'<span class="muted">ns</span>'}</td>`);
    tb.appendChild(tr);
  });
  t.appendChild(tb); s.appendChild(t); return s;
}
function rowDf(r){return r.df!=null?r.df:'';}

/* ---------- DIFFERENTIAL EXPRESSION ---------- */
let deState={cmp:null,ct:null};
function renderDE(sec){
  sec.innerHTML='';
  const d=ds();
  sec.appendChild(el('h2',null,'Differential expression — per cell type'));
  sec.appendChild(el('p','lead','Cell-type-resolved differential expression (Seurat Wilcoxon) for each pairwise contrast. Significant = FDR<0.05 and |log₂FC|≥0.25. Pick a contrast and cell type to explore the interactive volcano, top genes and the matched GO enrichment.'));

  const comps=Object.keys(d.de||{});
  if(!comps.length){sec.appendChild(el('div','empty','No per-cell-type DE results in this dataset.'));return;}

  // summary: significant genes per cell type, all comparisons
  const sumCard=card('Significant DE genes per cell type','Up- (▲ in first group) and down-regulated counts across all contrasts. Toggle a contrast in the legend.');
  plotInto(sumCard,'de-summary',Math.max(320,140+18*maxCT(d.de)));
  sec.appendChild(sumCard);

  // controls
  const ctrl=el('div','controls');
  const selCmp=el('select'); comps.forEach(c=>selCmp.appendChild(new Option(cmpLabel(c),c)));
  const selCt=el('select');
  ctrl.appendChild(el('label',null,'Contrast')); ctrl.appendChild(selCmp);
  ctrl.appendChild(el('label',null,'Cell type')); ctrl.appendChild(selCt);
  const dc=card('Volcano &amp; top differentially expressed genes');
  dc.appendChild(ctrl);
  const grid=el('div','grid2');
  const vWrap=el('div'); vWrap.id='de-volcano-wrap';
  const tWrap=el('div'); tWrap.id='de-table-wrap';
  grid.appendChild(vWrap); grid.appendChild(tWrap);
  dc.appendChild(grid);
  sec.appendChild(dc);

  const enrCard=card('Matched GO Biological-Process enrichment','Top over-represented terms for up- and down-regulated genes of the selected cell type / contrast (ORA, clusterProfiler).');
  enrCard.id='de-enr'; sec.appendChild(enrCard);

  function fillCts(){
    const cmp=selCmp.value; selCt.innerHTML='';
    Object.keys(d.de[cmp]).sort().forEach(ct=>selCt.appendChild(new Option(pretty(ct),ct)));
  }
  function update(){
    deState={cmp:selCmp.value,ct:selCt.value};
    drawVolcano('de-volcano-wrap',d.de[deState.cmp][deState.ct],deState.cmp,deState.ct,false);
    fillDeTable('de-table-wrap',d.de[deState.cmp][deState.ct]);
    fillDeEnr(d.enrichment.de,deState.cmp,deState.ct);
  }
  selCmp.onchange=()=>{fillCts();update();};
  selCt.onchange=update;
  fillCts(); update();
  drawDeSummary('de-summary',d.de);
}
function maxCT(de){let m=0;Object.values(de).forEach(b=>m=Math.max(m,Object.keys(b).length));return m*Object.keys(de).length;}
function drawDeSummary(id,de){
  const comps=Object.keys(de);
  // build per (comparison) traces of up and down stacked, grouped by cell type label "ct (cmp)"
  const rows=[];
  comps.forEach(cmp=>Object.keys(de[cmp]).forEach(ct=>rows.push({cmp,ct,...de[cmp][ct]})));
  rows.sort((a,b)=>(a.up+a.dn)-(b.up+b.dn));
  const labels=rows.map(r=>`${pretty(r.ct)}  ·  ${cmpLabel(r.cmp)}`);
  const tUp={x:rows.map(r=>r.up),y:labels,name:'Up in first group',orientation:'h',type:'bar',
    marker:{color:'#e0533d'},hovertemplate:'%{y}<br>Up: %{x}<extra></extra>'};
  const tDn={x:rows.map(r=>-r.dn),y:labels,name:'Down in first group',orientation:'h',type:'bar',
    marker:{color:'#3f7fd6'},customdata:rows.map(r=>r.dn),hovertemplate:'%{y}<br>Down: %{customdata}<extra></extra>'};
  Plotly.newPlot(id,[tDn,tUp],lay({barmode:'relative',
    xaxis:{title:'← down   significant genes (FDR<0.05, |log₂FC|≥0.25)   up →',zeroline:true,zerolinecolor:'#52608a'},
    yaxis:{automargin:true},margin:{l:230,r:14,t:10,b:45},legend:{orientation:'h',y:1.06}}),CFG);
}
function drawVolcano(wrapId,obj,cmp,ct,isPB){
  const wrap=$('#'+wrapId); wrap.innerHTML='';
  const h=el('h3',null,`Volcano — ${pretty(ct)}`); h.style.fontSize='14px';
  wrap.appendChild(h);
  wrap.appendChild(el('div','hint',`${cmpLabel(cmp)} · ${obj.up} up / ${obj.dn} down significant of ${obj.n_tested} tested`));
  const pdiv=el('div','plot'); pdiv.id=wrapId+'-p'; pdiv.style.minHeight='360px'; wrap.appendChild(pdiv);
  const pts=obj.volcano||[];
  const lfcT=isPB?1.0:0.25;
  const cat=p=>{const sig=p.padj!=null&&p.padj<0.05&&Math.abs(p.lfc)>=lfcT;
    if(!sig)return'ns'; return p.lfc>0?'up':'dn';};
  const groups={up:[],dn:[],ns:[]};
  pts.forEach(p=>{if(p.padj==null||p.padj<=0)return;groups[cat(p)].push(p);});
  const mk=(arr,color,name)=>({x:arr.map(p=>p.lfc),y:arr.map(p=>-Math.log10(p.padj)),
    mode:'markers',type:'scattergl',name,marker:{size:5,color,opacity:.6},
    text:arr.map(p=>`${p.gene}<br>log₂FC=${fmtN(p.lfc)} · FDR=${fmtP(p.padj)}`),
    hovertemplate:'%{text}<extra></extra>'});
  const traces=[mk(groups.ns,'#5a6680','ns'),
    mk(groups.dn,'#3f7fd6','Down'),mk(groups.up,'#e0533d','Up')];
  // label a few top genes
  const top=[...groups.up,...groups.dn].sort((a,b)=>a.padj-b.padj).slice(0,10);
  const ann=top.map(p=>({x:p.lfc,y:-Math.log10(p.padj),text:p.gene,showarrow:false,
    font:{size:9,color:'#e8edf6'},yshift:8}));
  Plotly.newPlot(pdiv.id,traces,lay({xaxis:{title:'log₂ fold-change',zeroline:true,zerolinecolor:'#52608a'},
    yaxis:{title:'−log₁₀ FDR'},shapes:[{type:'line',x0:lfcT,x1:lfcT,yref:'paper',y0:0,y1:1,line:{dash:'dot',color:'#52608a',width:1}},
    {type:'line',x0:-lfcT,x1:-lfcT,yref:'paper',y0:0,y1:1,line:{dash:'dot',color:'#52608a',width:1}}],
    annotations:ann,margin:{l:55,r:12,t:10,b:45},legend:{orientation:'h',y:1.08}}),CFG);
}
function fillDeTable(wrapId,obj){
  const wrap=$('#'+wrapId); wrap.innerHTML='';
  const h=el('h3',null,'Top genes'); h.style.fontSize='14px'; wrap.appendChild(h);
  wrap.appendChild(el('div','hint','Ranked by FDR. ▲ up / ▼ down in the first group of the contrast.'));
  const s=el('div','tbl-scroll'); const t=el('table','tbl');
  t.innerHTML='<thead><tr><th>Gene</th><th>Dir</th><th class="num">log₂FC</th><th class="num">FDR</th><th class="num">pct.1</th><th class="num">pct.2</th></tr></thead>';
  const tb=el('tbody');
  (obj.top||[]).forEach(g=>{
    const up=g.lfc>0;
    tb.appendChild(el('tr',null,`<td class="gene-name">${g.gene}</td>
      <td><span class="tag ${up?'up':'dn'}">${up?'▲':'▼'}</span></td>
      <td class="num ${up?'dirpos':'dirneg'}">${fmtN(g.lfc)}</td>
      <td class="num">${fmtP(g.padj)}</td>
      <td class="num">${g.pct1!=null?g.pct1.toFixed(2):'—'}</td>
      <td class="num">${g.pct2!=null?g.pct2.toFixed(2):'—'}</td>`));
  });
  if(!(obj.top||[]).length) tb.appendChild(el('tr',null,'<td colspan="6" class="empty">No significant genes.</td>'));
  t.appendChild(tb); s.appendChild(t); wrap.appendChild(s);
}
function fillDeEnr(enrDe,cmp,ct){
  const host=$('#de-enr');
  // keep title/hint (first two children), drop rest
  while(host.children.length>2) host.removeChild(host.lastChild);
  const d=(enrDe&&enrDe[cmp]&&enrDe[cmp][ct])||null;
  const grid=el('div','grid2');
  grid.appendChild(termsBlock('Enriched in UP genes',d&&d.GO_BP_UP,'#e0533d'));
  grid.appendChild(termsBlock('Enriched in DOWN genes',d&&d.GO_BP_DOWN,'#3f7fd6'));
  host.appendChild(grid);
  // dotplot figures
  const figs=el('div','gallery');
  [['UP','GO-BP enrichment (up)'],['DOWN','GO-BP enrichment (down)']].forEach(([dir,cap])=>{
    const t=thumb(fig(`enr_dot_${cmp}_${ct}_${dir}`),cap);if(t)figs.appendChild(t);
  });
  if(figs.children.length) host.appendChild(figs);
}
function termsBlock(title,terms,color){
  const b=el('div');
  b.appendChild(el('h4',null,title));
  if(!terms||!terms.length){b.appendChild(el('div','empty','No enriched terms.'));return b;}
  const max=Math.max(...terms.map(t=>-Math.log10(t.padj||1)));
  const ul=el('ul','terms');
  terms.forEach(t=>{
    const w=Math.max(8,80*(-Math.log10(t.padj||1)/(max||1)));
    const li=el('li');
    li.innerHTML=`<span class="bar" style="width:${w}px;background:${color}"></span>
      <span class="desc">${t.desc}<br><span class="q">${t.id} · ${t.count} genes · FDR=${fmtP(t.padj)}</span></span>`;
    ul.appendChild(li);
  });
  b.appendChild(ul); return b;
}

/* ---------- PSEUDOBULK ---------- */
function renderPseudobulk(sec){
  sec.innerHTML='';
  const d=ds();
  sec.appendChild(el('h2',null,'Pseudo-bulk differential expression'));
  sec.appendChild(el('p','lead','Replicate-level pseudo-bulk DESeq2 contrasts aggregate counts per sample, giving a conservative, biologically-replicated view of expression change. Significant = FDR<0.05 and |log₂FC|≥1.'));
  const comps=Object.keys(d.pseudobulk||{});
  if(!comps.length){sec.appendChild(el('div','empty','No pseudo-bulk results in this dataset.'));return;}

  const ctrl=el('div','controls');
  const sel=el('select'); comps.forEach(c=>sel.appendChild(new Option(cmpLabel(c),c)));
  ctrl.appendChild(el('label',null,'Contrast')); ctrl.appendChild(sel);
  const pc=card('Pseudo-bulk volcano &amp; top genes');
  pc.appendChild(ctrl);
  const grid=el('div','grid2');
  const v=el('div'); v.id='pb-volcano-wrap';
  const t=el('div'); t.id='pb-table-wrap';
  grid.appendChild(v); grid.appendChild(t); pc.appendChild(grid);
  sec.appendChild(pc);

  const enr=card('Pseudo-bulk GO Biological-Process enrichment','Top over-represented terms among up/down pseudo-bulk genes.');
  enr.id='pb-enr'; sec.appendChild(enr);

  function update(){
    const cmp=sel.value, o=d.pseudobulk[cmp];
    drawVolcano('pb-volcano-wrap',o,cmp,'pseudo-bulk',true);
    fillPbTable('pb-table-wrap',o);
    // enrichment
    while(enr.children.length>2)enr.removeChild(enr.lastChild);
    const e=d.enrichment.pseudo_bulk&&d.enrichment.pseudo_bulk[cmp];
    const g=el('div','grid2');
    g.appendChild(termsBlock('Enriched in UP genes',e&&e.GO_BP_UP,'#e0533d'));
    g.appendChild(termsBlock('Enriched in DOWN genes',e&&e.GO_BP_DOWN,'#3f7fd6'));
    enr.appendChild(g);
    const figs=el('div','gallery');
    [['UP','GO-BP (up)'],['DOWN','GO-BP (down)']].forEach(([dir,cap])=>{const th=thumb(fig(`enr_pb_dot_${cmp}_${dir}`),cap);if(th)figs.appendChild(th);});
    const vol=thumb(fig('pb_volcano_'+cmp),'Pipeline volcano (DESeq2)');if(vol)figs.appendChild(vol);
    if(figs.children.length)enr.appendChild(figs);
  }
  sel.onchange=update; update();
}
function fillPbTable(wrapId,o){
  const wrap=$('#'+wrapId); wrap.innerHTML='';
  const h=el('h3',null,'Top genes'); h.style.fontSize='14px'; wrap.appendChild(h);
  wrap.appendChild(el('div','hint',`${o.up} up / ${o.dn} down significant of ${o.n_tested} tested. Ranked by FDR.`));
  const s=el('div','tbl-scroll'); const t=el('table','tbl');
  t.innerHTML='<thead><tr><th>Gene</th><th>Dir</th><th class="num">log₂FC</th><th class="num">FDR</th><th class="num">baseMean</th></tr></thead>';
  const tb=el('tbody');
  (o.top||[]).forEach(g=>{const up=g.lfc>0;
    tb.appendChild(el('tr',null,`<td class="gene-name">${g.gene}</td>
      <td><span class="tag ${up?'up':'dn'}">${up?'▲':'▼'}</span></td>
      <td class="num ${up?'dirpos':'dirneg'}">${fmtN(g.lfc)}</td>
      <td class="num">${fmtP(g.padj)}</td>
      <td class="num">${g.baseMean!=null?g.baseMean.toFixed(0):'—'}</td>`));
  });
  if(!(o.top||[]).length)tb.appendChild(el('tr',null,'<td colspan="5" class="empty">No significant genes.</td>'));
  t.appendChild(tb);s.appendChild(t);wrap.appendChild(s);
}

/* ---------- ENRICHMENT ---------- */
function renderEnrichment(sec){
  sec.innerHTML='';
  const d=ds();
  sec.appendChild(el('h2',null,'Functional enrichment (GO over-representation)'));
  sec.appendChild(el('p','lead','Over-representation analysis (clusterProfiler, org.Mm.eg.db) of differentially-expressed gene sets against Gene Ontology Biological Process, Cellular Component and Molecular Function. Browse any contrast / cell type, split by direction of change.'));

  const enrDe=d.enrichment.de||{};
  const comps=Object.keys(enrDe);
  if(!comps.length){sec.appendChild(el('div','empty','No enrichment results in this dataset.'));return;}

  const ctrl=el('div','controls');
  const selCmp=el('select'); comps.forEach(c=>selCmp.appendChild(new Option(cmpLabel(c),c)));
  const selCt=el('select');
  const selOnto=el('select'); [['GO_BP','GO: Biological Process'],['GO_CC','GO: Cellular Component'],['GO_MF','GO: Molecular Function']].forEach(([v,l])=>selOnto.appendChild(new Option(l,v)));
  ctrl.appendChild(el('label',null,'Contrast'));ctrl.appendChild(selCmp);
  ctrl.appendChild(el('label',null,'Cell type'));ctrl.appendChild(selCt);
  ctrl.appendChild(el('label',null,'Ontology'));ctrl.appendChild(selOnto);
  const c=card('Enriched terms by direction');
  c.appendChild(ctrl);
  const host=el('div'); host.id='enr-host'; c.appendChild(host);
  sec.appendChild(c);

  function fillCts(){selCt.innerHTML='';Object.keys(enrDe[selCmp.value]).sort().forEach(ct=>selCt.appendChild(new Option(pretty(ct),ct)));}
  function update(){
    const cmp=selCmp.value,ct=selCt.value,onto=selOnto.value;
    const data=enrDe[cmp][ct];
    host.innerHTML='';
    const g=el('div','grid2');
    g.appendChild(termsBlock(`Up — ${onto.replace('GO_','GO ')}`,data[onto+'_UP'],'#e0533d'));
    g.appendChild(termsBlock(`Down — ${onto.replace('GO_','GO ')}`,data[onto+'_DOWN'],'#3f7fd6'));
    host.appendChild(g);
    if(onto==='GO_BP'){
      const figs=el('div','gallery');
      [['UP','GO-BP dot-plot (up)'],['DOWN','GO-BP dot-plot (down)']].forEach(([dir,cap])=>{const t=thumb(fig(`enr_dot_${cmp}_${ct}_${dir}`),cap);if(t)figs.appendChild(t);});
      if(figs.children.length)host.appendChild(figs);
    }
  }
  selCmp.onchange=()=>{fillCts();update();};
  selCt.onchange=update; selOnto.onchange=update;
  fillCts(); update();
}

/* ---------- SPLICING ---------- */
function renderSplicing(sec){
  sec.innerHTML='';
  const d=ds();
  sec.appendChild(el('h2',null,'Splicing ratios'));
  sec.appendChild(el('p','lead','Per-cell splicing ratio = unspliced ⁄ (spliced + 1), derived from spliced/unspliced count layers. Elevated ratios indicate accumulation of unspliced transcript — the expected signature of spliceosome disruption. Compared globally and per cell type across conditions.'));

  sec.appendChild(el('div','callout',`<b>Both mutants elevate the splicing ratio.</b> Control = ${fmtN(d.splicing_global.Control)}, <i>Eftud2</i> = ${fmtN(d.splicing_global.Eftud2)}, <i>Sf3b4</i> = ${fmtN(d.splicing_global.Sf3b4)} (n-weighted global means).`));

  const g=el('div','grid2');
  const gc=card('Global splicing ratio by condition','n-weighted mean across cell types.');
  plotInto(gc,'spl-global',360); g.appendChild(gc);
  const cc=card('Splicing ratio by cell type &amp; condition','Mean ± SD per cell type. Hover for n.');
  plotInto(cc,'spl-ct',360); g.appendChild(cc);
  sec.appendChild(g);

  const figs=card('Reference splicing figures','Pipeline-generated distributions.');
  figs.appendChild(gallery([['splicing_by_condition','Splicing ratio distribution by condition'],
    ['splicing_by_celltype','Splicing ratio by cell type']]));
  sec.appendChild(figs);

  drawSplicingGlobal('spl-global');
  drawSplicingByCT('spl-ct');
}
function drawSplicingByCT(id){
  const s=ds().splicing;
  const cts=[...new Set(s.map(r=>r.cell_type))];
  const traces=COND.map(cn=>{
    const rows=cts.map(ct=>s.find(r=>r.cell_type===ct&&r.condition===cn));
    return {x:cts.map(pretty),y:rows.map(r=>r?r.mean:null),name:cn,type:'bar',marker:{color:COL[cn]},
      error_y:{type:'data',array:rows.map(r=>r?r.sd:null),visible:true,color:'#6b7794',thickness:1,width:2},
      customdata:rows.map(r=>r?r.n:0),
      hovertemplate:`${cn} · %{x}<br>ratio=%{y:.2f}±%{error_y.array:.2f}<br>n=%{customdata}<extra></extra>`};
  });
  Plotly.newPlot(id,traces,lay({barmode:'group',yaxis:{title:'unspliced/(spliced+1)'},
    xaxis:{tickangle:-30},margin:{l:60,r:10,t:10,b:110},legend:{orientation:'h',y:1.1}}),CFG);
}

/* ---------- VELOCITY ---------- */
function renderVelocity(sec){
  sec.innerHTML='';
  const d=ds();
  sec.appendChild(el('h2',null,'RNA velocity dynamics'));
  sec.appendChild(el('p','lead','Summary scVelo metrics per condition — velocity speed (magnitude), coherence (directional consistency) and mean latent time — capturing how transcriptional dynamics are reshaped in the mutants.'));
  const v=d.velocity;
  if(!v||!v.length){sec.appendChild(el('div','empty','No velocity summary available.'));return;}
  const g=el('div','grid2');
  const bc=card('Velocity metrics by condition','Grouped comparison of the three summary statistics.');
  plotInto(bc,'vel-bars',380); g.appendChild(bc);
  const rc=card('Condition velocity profile','Radar of normalised metrics.');
  plotInto(rc,'vel-radar',380); g.appendChild(rc);
  sec.appendChild(g);

  const figs=card('Reference velocity figures','Pipeline outputs.');
  figs.appendChild(gallery([['08_velocity_summary_bars','Velocity summary bars'],
    ['06_velocity_metrics_tile_heatmap','Velocity metrics heatmap'],
    ['07_delta_velocity_scatter','Δ velocity scatter'],
    ['10_shared_genes_deltavelocity_heatmap','Shared Δ-velocity genes'],
    ['09_sig_genes_per_celltype','Significant genes per cell type']]));
  sec.appendChild(figs);

  // grouped bars
  const metrics=[['speed','Speed'],['coherence','Coherence'],['latent_time','Latent time']];
  const traces=v.map(r=>({x:metrics.map(m=>m[1]),y:metrics.map(m=>r[m[0]]),name:r.condition,type:'bar',
    marker:{color:COL[r.condition]},hovertemplate:`${r.condition} · %{x}: %{y:.3f}<extra></extra>`}));
  Plotly.newPlot('vel-bars',traces,lay({barmode:'group',yaxis:{title:'value'},margin:{l:55,r:10,t:10,b:40},legend:{orientation:'h',y:1.1}}),CFG);
  // radar normalised 0..1 per metric
  const norm=m=>{const vals=v.map(r=>r[m]);const mn=Math.min(...vals),mx=Math.max(...vals);return r=>mx>mn?(r[m]-mn)/(mx-mn):0.5;};
  const ns=metrics.map(m=>norm(m[0]));
  const rt=v.map(r=>({type:'scatterpolar',r:metrics.map((m,i)=>ns[i](r)),theta:metrics.map(m=>m[1]),
    fill:'toself',name:r.condition,line:{color:COL[r.condition]},
    marker:{color:COL[r.condition]},opacity:.65}));
  Plotly.newPlot('vel-radar',rt,lay({polar:{bgcolor:'rgba(0,0,0,0)',
    radialaxis:{visible:true,range:[0,1],gridcolor:'#243049',tickfont:{size:9}},
    angularaxis:{gridcolor:'#243049'}},margin:{l:40,r:40,t:30,b:30},legend:{orientation:'h',y:1.12}}),CFG);
}

/* ---------- METHODS ---------- */
function renderMethods(sec){
  sec.innerHTML='';
  sec.appendChild(el('h2',null,'Methods &amp; data provenance'));
  sec.appendChild(el('p','lead','Summary of the computational workflow. Full, runnable R pipelines are versioned in the repository.'));
  const c=card('Pipeline summary');
  c.innerHTML+=`<div style="font-size:13.5px;line-height:1.75;color:#c7d1e6">
    <p><b>Design.</b> E9.5 mouse embryos, three conditions — wild-type <b>Control</b> and two craniofacial spliceosomopathy mutants, <i>Eftud2</i> and <i>Sf3b4</i> — each pooling two biological replicates (10x Genomics; spliced/unspliced quantified for RNA velocity).</p>
    <p><b>Processing.</b> Per-sample loading of spliced & unspliced layers, QC, normalisation, integration, dimensionality reduction (PCA → UMAP/t-SNE), graph clustering and cell-type annotation (Seurat v3 assay model). Two atlases are analysed: the full embryo (<b>all cell types</b>) and a re-clustered <b>neural-crest</b> compartment.</p>
    <p><b>Composition.</b> Per-condition cell-type fractions; significance via pooled χ², per-cell-type Fisher exact tests (odds ratios, BH-adjusted) and a replicate-aware centred-log-ratio linear mixed model (omnibus LRT + pairwise contrasts).</p>
    <p><b>Differential expression.</b> Cell-type-resolved Wilcoxon tests (Seurat <span class="mono">FindMarkers</span>) for the three pairwise contrasts; significance at FDR&lt;0.05, |log₂FC|≥0.25. Replicate-level <b>pseudo-bulk</b> contrasts with DESeq2 (FDR&lt;0.05, |log₂FC|≥1).</p>
    <p><b>Functional enrichment.</b> Over-representation analysis (clusterProfiler, org.Mm.eg.db) against GO BP/CC/MF, split by direction; dot-plots and enrichment maps per gene set.</p>
    <p><b>Splicing.</b> Per-cell splicing ratio = unspliced ⁄ (spliced + 1); summarised globally and per cell type per condition.</p>
    <p><b>RNA velocity.</b> scVelo-based velocity (via velociraptor/reticulate); per-condition speed, coherence and latent-time summaries plus differential-velocity genes.</p>
    <p><b>Software.</b> Seurat, SingleCellExperiment, scran/scater, velociraptor/scVelo, DESeq2, limma, clusterProfiler, enrichplot, EnhancedVolcano, ggplot2.</p>
  </div>`;
  sec.appendChild(c);
  const c2=card('Conditions, contrasts &amp; thresholds');
  c2.innerHTML+=`<div class="pillrow">
    <span class="badge">Control = baseline</span><span class="badge">Eftud2 vs Control</span>
    <span class="badge">Sf3b4 vs Control</span><span class="badge">Eftud2 vs Sf3b4</span>
    <span class="badge">DE: FDR&lt;0.05 · |log₂FC|≥0.25</span><span class="badge">Pseudo-bulk: FDR&lt;0.05 · |log₂FC|≥1</span></div>
    <div class="muted" style="font-size:13px">Scripts: <span class="mono">process_ALL_CELLS.R</span>, <span class="mono">process_NEURAL_CREST_CELLS.R</span>. See the repository README for the full description, citation and license.</div>`;
  sec.appendChild(c2);
  const c3=card('Citation &amp; license');
  c3.innerHTML+=`<div style="font-size:13.5px;color:#c7d1e6">
    <p>Unpublished research data — <b>restricted use</b>. Permitted users: Prof. Loydie A. Majewska, Dr. Joshua Moore, Dr. Marie-Claude Beauchamp and Dr. Eric Bareke (McGill University). Please cite as described in <span class="mono">CITATION.cff</span> and do not redistribute figures or data without permission.</p></div>`;
  sec.appendChild(c3);
}

/* ============================================================================
 * ROUTER
 * ==========================================================================*/
const RENDER={overview:renderOverview,composition:renderComposition,de:renderDE,
  pseudobulk:renderPseudobulk,enrichment:renderEnrichment,splicing:renderSplicing,
  velocity:renderVelocity,methods:renderMethods};
function renderAll(){
  document.querySelectorAll('.section').forEach(s=>{
    const f=RENDER[s.dataset.sec]; if(f) f(s);
  });
  setTimeout(()=>window.dispatchEvent(new Event('resize')),80);
}
function setDataset(ds_){
  DS=ds_;
  document.querySelectorAll('#dsToggle button').forEach(b=>b.classList.toggle('active',b.dataset.ds===ds_));
  renderAll();
}
$('#dsToggle').addEventListener('click',e=>{const b=e.target.closest('button');if(b)setDataset(b.dataset.ds);});

buildNav();
renderAll();
