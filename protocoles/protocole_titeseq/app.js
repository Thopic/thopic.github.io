/* ---------------- constants & helpers ---------------- */
const NA = 6.022e23;
const LADDER = [200,400,600,1000,1400,2000,3000,6000,10000,20000,50000,100000,200000,500000];
const $ = id => document.getElementById(id);
const el = (t,c,h) => { const n=document.createElement(t); if(c)n.className=c; if(h!=null)n.innerHTML=h; return n; };

function sci(x,d){ if(x===0)return '0'; if(!isFinite(x))return '—';
  const e=Math.floor(Math.log10(Math.abs(x)));
  return (x/Math.pow(10,e)).toFixed(d==null?2:d)+'e'+(e<0?'':'+')+e; }
function vol(u){ if(u==null||!isFinite(u))return '—';
  if(u>=1000) return (u/1000).toFixed(u>=10000?0:1).replace(/[.,]0$/,'').replace('.',',')+' mL';
  if(u>=10) return u.toFixed(0)+' µL';
  if(u>=1) return u.toFixed(2)+' µL';
  return u.toFixed(3)+' µL'; }
function cells(n){ if(!isFinite(n))return '—';
  if(n>=1e6) return (n/1e6).toFixed(n>=1e7?1:2)+'e6';
  if(n>=1e3) return (n/1e3).toFixed(0)+'e3';
  return n.toFixed(0); }
function vessel(v){ if(v<=1500)return 'tube 1,5 mL'; if(v<=15000)return 'conique 15 mL';
  if(v<=50000)return 'conique 50 mL'; if(v<=500000)return 'conique 500 mL'; return 'à répartir sur plusieurs récipients'; }

/* stock de ligand : molarité en µM, quelle que soit l'unité fournie */
function stockUM(){ return S.stockMode==='mass' ? 1000*(+S.stockVal)/(+S.mw||1) : (+S.stockVal); }
function stockMgMl(){ return stockUM()*(+S.mw)/1000; }
function massUG(uL){ return uL*stockUM()*(+S.mw)/1000; }
function mass(ug){ return ug>=1000 ? (ug/1000).toFixed(2)+' mg' : ug.toFixed(ug<10?2:1)+' µg'; }
function worstKd(){ return Math.pow(10,-(+S.worstLogKd)); }
function LIG(){ return S.ligName||'ligand'; }
function DISP(){ return S.dispName||'protéine'; }
function cap(t){ return t.charAt(0).toUpperCase()+t.slice(1); }

/* 1:1 binding — complex concentration */
function complexConc(A,R,Kd){ const b=A+R+Kd; const d=b*b-4*A*R;
  return (b-Math.sqrt(Math.max(d,0)))/2; }

/* ---------------- state ---------------- */
const DEFAULT_VOLS = {'7':200,'7.5':600,'8':1400,'8.5':6000,'9':20000,'9.5':50000,
  '10':200000,'10.5':500000,'11':500000,'11.5':500000,'12':500000,'12.5':500000};
const DEFAULT_CELLS = {'12.5':500};

const DEF = {
  anchor:'', indH:16, primH:16, statH:22,
  ligName:'ligand', dispName:'protéine', libName:'LIB',
  mw:85, stockMode:'molar', stockVal:3.5, stockHave:50, cellsPerOD:1e7, recPerCell:5e4, worstLogKd:11, maxDepl:10, minPip:1,
  diversity:65536, cellVol:700, cellVolCtl:50, sortTarget:1e6,
  lossLabel:50, lossWaste:30, lossRecov:30, nPEbins:4, nFITCbins:8,
  zeroVol:500, minHeadroom:50, logStart:7, logEnd:12, logStep:1, maxVessel:500000, cloneNames:'CTL1,CTL2',
  cloneConc:['8','10','12'], cloneVols:{'8':200,'10':10000,'12':50000},
  vols:Object.assign({},DEFAULT_VOLS), cellVols:Object.assign({},DEFAULT_CELLS),
  satOD:3.9, indOD:0.67, indVol:5,
  odPrep:{}, libResus:12000, cloneResus:1000,
  secVol:250, secOver:10, dilPE:100, dilMyc:50, dilFITC:50,
  seed:1, sortTime:{}, ctlEvents:30000, doPlating:true, targetCol:500, plateSamples:'9.5,10', outVol:4, targetCol2:500, qcConc:'9.5',
  sorted:{}, binCells:{}, plateFlags:{}, odLog:{},
  counts:[{s:'pré-tri e9.5',y:'',d:'',e:''},{s:'pré-tri e10',y:'',d:'',e:''}]
};
let S = JSON.parse(JSON.stringify(DEF));

const STORE_KEY='titeseq-rbd-fr-v1';
function save(){ try{ localStorage.setItem(STORE_KEY, JSON.stringify(S)); }catch(e){} }
function load(){ try{ const r=localStorage.getItem(STORE_KEY);
  if(r) S=Object.assign(JSON.parse(JSON.stringify(DEF)), JSON.parse(r)); }catch(e){} }

/* ---------------- derived model ---------------- */
function ladder(){
  const out=[]; const a=+S.logStart, b=+S.logEnd, st=Math.abs(+S.logStep)||0.5;
  for(let x=a; x<=b+1e-9; x+=st) out.push(Math.round(x*1000)/1000);
  return out;
}
function key(l){ return String(l); }
/* reproducible shuffle: the seed is recorded so the day's order can be reconstructed */
function rng(seed){ let a=(seed>>>0)||1;
  return ()=>{ a=a+0x6D2B79F5|0; let t=Math.imul(a^a>>>15,1|a);
    t=t+Math.imul(t^t>>>7,61|t)^t; return ((t^t>>>14)>>>0)/4294967296; }; }
function runOrder(){
  const L=ladder().slice(), r=rng(+S.seed||1);
  for(let i=L.length-1;i>0;i--){ const j=Math.floor(r()*(i+1)); const x=L[i]; L[i]=L[j]; L[j]=x; }
  return L;
}
function clones(){ return String(S.cloneNames).split(',').map(s=>s.trim()).filter(Boolean); }

function rowFor(l){
  const k=key(l);
  const A=Math.pow(10,-l);
  const cv=S.cellVols[k]!=null ? +S.cellVols[k] : +S.cellVol;
  const tv=S.vols[k]!=null ? +S.vols[k] : autoVolume(l,cv);
  const nCells=cv/1000*(+S.cellsPerOD);
  const R=nCells*(+S.recPerCell)/NA/(tv*1e-6);
  const C=complexConc(A,R,worstKd());
  return {l,k,A,cv,tv,nCells,R,ratio:A/R,free:A-C,depl:100*C/A,
    isClone:S.cloneConc.includes(k)};
}
function autoVolume(l,cv){
  const A=Math.pow(10,-l), cap=+S.maxVessel;
  const cand=LADDER.filter(v=>v<=cap);
  for(const v of cand){
    const R=(cv/1000*(+S.cellsPerOD))*(+S.recPerCell)/NA/(v*1e-6);
    if(100*complexConc(A,R,worstKd())/A <= +S.maxDepl) return v;
  }
  return cap;
}
/* total volume of primary solution to make at a concentration */
/* cells go straight into the solution when there is room for them;
   otherwise the sample is pelleted and taken up in the full volume */
function addMode(total, cv){
  const head=total-cv;
  if(head>=(+S.minHeadroom)) return {direct:true, label:'ajout direct de '+cv+' µL'};
  return {direct:false, label:'culot'};
}

/* the library reaction only; control solutions are prepared separately */
function primaryTotal(r){ return {lib:r.tv, total:r.tv}; }

/* every control reaction, grouped by the solution it shares */
function controlRxns(){
  const out=[], L=ladder(), cl=clones(), LB=S.libName||'LIB';
  if(L.length) out.push({kind:'comp', label:LB+'_'+LIG()+'PE', strains:['banque'],
    l:L[0], n:1, unit:+S.cellVolCtl, cv:+S.cellVolCtl,
    why:'compensation PE'});
  S.cloneConc.filter(k=>L.map(key).includes(k)).forEach(k=>{
    out.push({kind:'clone', label:cl.join(', '), strains:cl, l:parseFloat(k), n:cl.length,
      unit:+S.cloneVols[k]||0, cv:+S.cellVolCtl, why:'contrôle clonal'});
  });
  return out.filter(r=>r.n>0&&r.unit>0);
}
function ctrlRow(c){
  const A=Math.pow(10,-c.l), total=c.n*c.unit;
  const R=(c.cv/1000*(+S.cellsPerOD))*(+S.recPerCell)/NA/(c.unit*1e-6);
  const C=complexConc(A,R,worstKd());
  return {A, total, R, depl:100*C/A};
}
/* ACE2 source cascade */
function aceSource(totalUL, A){
  const stockM=stockUM()*1e-6;
  const need = totalUL*A/stockM;           // µL of neat stock
  const minP = +S.minPip;
  let f=1;
  if(need<minP) f=10;
  if(need*10<minP) f=100;
  return {factor:f, addUL:need*f, neatEquiv:need,
    label: f===1?'pur':'1:'+f};
}
function sampleList(){
  const out=[]; const L=ladder(); const top=L[0];
  const cvCtl=+S.cellVolCtl, cv=+S.cellVol;
  const LB=S.libName||'LIB', zv=+S.zeroVol;
  out.push({n:LB+'_unlabeled',strain:'banque',prim:'—',sec:'—',cv:cv,tv:zv,why:'fenêtre FSC/SSC, compensation'});
  out.push({n:LB+'_mycFITC',strain:'banque',prim:'—',sec:'myc-FITC',cv:cv*2,tv:zv,why:'tri d’expression, compensation'});
  out.push({n:LB+'_'+LIG()+'PE',strain:'banque',prim:'e'+top,sec:'strep-PE',cv:cvCtl,tv:cvCtl,why:'compensation PE'});
  out.push({n:LB+'_FITC',strain:'banque',prim:'—',sec:'FITC',cv:cvCtl,tv:zv,why:'fenêtre de bruit de fond FITC'});
  out.push({n:LB+'_0',strain:'banque',prim:'—',sec:'strep-PE + myc-FITC',cv:cv,tv:zv,why:'fenêtre de bruit de fond PE'});
  ladder().forEach(l=>{ const r=rowFor(l);
    out.push({n:LB+'_'+l,strain:'banque',prim:'e'+l,sec:'strep-PE + myc-FITC',cv:r.cv,tv:r.tv,
      why:'tri de liaison',conc:l}); });
  clones().forEach(c=>{
    out.push({n:c+'_noACE',strain:c,prim:'—',sec:'strep-PE + myc-FITC',cv:cvCtl,tv:zv,why:'contrôle clonal'});
    S.cloneConc.filter(k=>ladder().map(key).includes(k)).forEach(k=>{
      out.push({n:c+'_'+k,strain:c,prim:'e'+k,sec:'strep-PE + myc-FITC',cv:cvCtl,tv:+S.cloneVols[k]||0,why:'contrôle clonal'}); });
  });
  return out;
}
/* the order samples are run on the instrument: controls first, then the shuffled ladder */
function orderedSamples(){
  const sl=sampleList(), ord=runOrder();
  const rank=s=>{
    if(s.conc!=null) return 100+ord.indexOf(s.conc);
    if(s.strain!=='banque') return 1000+sl.indexOf(s);
    return sl.indexOf(s);
  };
  return sl.map((s,i)=>({s,i,r:rank(s)})).sort((a,b)=>a.r-b.r||a.i-b.i).map(x=>x.s);
}
function budget(){
  const labeled=(+S.cellVol)/1000*(+S.cellsPerOD);
  const survive=labeled*(1-(+S.lossLabel)/100);
  const sorted=Math.min(+S.sortTarget, survive);
  const inBins=sorted*(1-(+S.lossWaste)/100);
  const recovered=inBins*(1-(+S.lossRecov)/100);
  const nBind=ladder().length;
  return {labeled,survive,sorted,inBins,recovered,cover:recovered/(+S.diversity),
    totalSorted:sorted*nBind, nBind};
}

/* ---------------- rendering ---------------- */
function fmtDate(d){ if(!d||isNaN(d))return '—';
  return d.toLocaleString('fr-FR',{weekday:'short',month:'numeric',day:'numeric',
    hour:'numeric',minute:'2-digit'}); }

function renderClock(){
  const t=$('timeline'); t.innerHTML='';
  const base=S.anchor?new Date(S.anchor):null;
  const h=x=>base?new Date(base.getTime()+x*3600e3):null;
  const steps=[
    ['Décongélation de la banque', -(+S.indH)-(+S.statH), 'début de la phase stationnaire'],
    ['Induction en SGDCAA', -(+S.indH), (+S.statH)+' h après décongélation'],
    ['Marquage primaire', 0, (+S.indH)+' h après induction'],
    ['Marquage secondaire', +S.primH, 'Après ' + (S.primH)+' h de primaire'],
    ['Début du tri', (+S.primH)+3, 'Après 45 min de secondaire + lavages'],
    ['Récupération et étalement', (+S.primH)+7, 'après ~4 h sur l’appareil'],
    ['Congélation à DO 0,9–2', (+S.primH)+26, 'vérifier toutes les 2 h dès le matin']
  ];
  steps.forEach(([what,off,gap])=>{
    const d=el('div','tl');
    d.appendChild(el('div','when',fmtDate(h(off))));
    d.appendChild(el('div','what',what));
    d.appendChild(el('div','gap',gap));
    t.appendChild(d);
  });
  const tags={'t-day1':-(+S.indH)-(+S.statH),'t-day2':-(+S.indH),'t-day3':0,
    't-day4':+S.primH,'t-sort':(+S.primH)+3,'t-recov':(+S.primH)+10,'t-day5':(+S.primH)+26};
  Object.entries(tags).forEach(([id,off])=>{ const n=$(id); if(n) n.textContent=base?fmtDate(h(off)):'définir une heure de départ'; });
}

function renderBudget(){
  const b=budget();
  const cls=b.cover<10?'bad':(b.cover<20?'warn':'good');
  $('budget').innerHTML=
   '<div class="grid" style="margin-top:14px">'+
   kv('Cellules marquées par échantillon',cells(b.labeled))+
   kv('Survivant au marquage',cells(b.survive))+
   kv('Triées par échantillon',cells(b.sorted))+
   kv('Arrivant dans les bins',cells(b.inBins))+
   kv('Récupérées par échantillon','<span class="'+cls+'">'+cells(b.recovered)+'</span>')+
   kv('Couverture de la banque','<span class="'+cls+'">'+b.cover.toFixed(1)+'× la diversité</span>')+
   kv('Cellules totales sur '+b.nBind+' concentrations',cells(b.totalSorted))+
   '</div>'+
   '';
}
function kv(k,v){ return '<div class="kv"><span>'+k+'</span><span class="v">'+v+'</span></div>'; }

function renderCurve(){
  const t=$('curveTable');
  const L=LIG(), D=DISP();
  let h='<thead><tr><th class="l">Échantillon</th><th>−log['+L+']</th><th>['+L+'] (M)</th>'+
    '<th>Cellules (µL)</th><th>Vol. marquage (µL)</th><th class="l">Récipient</th>'+
    '<th>['+D+'] (M)</th><th>'+L+':'+D+'</th><th>'+cap(L)+' libre (M)</th><th>Déplétion</th>'+
    '</tr></thead><tbody>';
  ladder().forEach(l=>{
    const r=rowFor(l);
    const dcls=r.depl>(+S.maxDepl)?'badcell bad':(r.depl>(+S.maxDepl)*0.8?'warncell warn':'');
    h+='<tr><td class="l">'+(S.libName||'LIB')+'_'+l+'</td><td>'+l.toFixed(2)+'</td><td>'+sci(r.A)+'</td>'+
      '<td><input class="mini" type="number" step="50" data-cv="'+r.k+'" value="'+r.cv+'"></td>'+
      '<td><input class="mini wide" type="number" step="100" data-tv="'+r.k+'" value="'+r.tv+'"></td>'+
      '<td class="l">'+vessel(r.tv)+'</td><td>'+sci(r.R)+'</td><td>'+r.ratio.toFixed(1)+'</td>'+
      '<td>'+sci(r.free)+'</td><td class="'+dcls+'">'+r.depl.toFixed(1)+'%</td></tr>';
  });
  const rx=controlRxns();
  if(rx.length){
    h+='<tr class="sub"><td class="l" colspan="10">Contrôles clonaux — réactions séparées</td></tr>';
    rx.forEach(c=>{
      const cr=ctrlRow(c);
      const dcls=cr.depl>(+S.maxDepl)?'badcell bad':(cr.depl>(+S.maxDepl)*0.8?'warncell warn':'');
      h+='<tr><td class="l">'+c.label+' <span class="muted">· '+c.why+'</span></td>'+
        '<td>'+c.l.toFixed(2)+'</td><td>'+sci(cr.A)+'</td><td>'+c.cv+'</td>'+
        '<td>'+(c.kind==='clone'
          ? '<input class="mini wide" type="number" step="50" data-clv="'+key(c.l)+'" value="'+c.unit+'">'
          : c.unit.toFixed(0))+'</td>'+
        '<td class="l">'+vessel(cr.total)+(c.n>1?' <span class="muted">× '+c.n+'</span>':'')+'</td>'+
        '<td>'+sci(cr.R)+'</td><td>'+(cr.A/cr.R).toFixed(1)+'</td>'+
        '<td>'+sci(cr.A-complexConc(cr.A,cr.R,worstKd()))+'</td>'+
        '<td class="'+dcls+'">'+cr.depl.toFixed(1)+'%</td></tr>';
    });
  }
  h+='</tbody>';
  t.innerHTML=h;
  $('curveNotes').innerHTML='';
}

/* ligand drawn by every reaction on the page, library and controls alike */
function stockSummary(){
  let neat=0,d10=0,d100=0;
  const take=(v,A)=>{ const s=aceSource(v,A); neat+=s.neatEquiv;
    if(s.factor===10)d10+=s.addUL; if(s.factor===100)d100+=s.addUL; };
  ladder().forEach(l=>{ const r=rowFor(l); take(primaryTotal(r).total, r.A); });
  controlRxns().forEach(c=>{ take(ctrlRow(c).total, ctrlRow(c).A); });
  return {neat,d10,d100};
}
function renderPrimary(){
  const t=$('primeTable');
  let h='<thead><tr><th>−log</th><th>['+LIG()+'] (M)</th><th>Volume final (µL)</th>'+
    '<th class="l">Récipient</th><th class="l">Source</th><th>Ajouter (µL)</th>'+
    '<th>PBSA-BL (µL)</th><th>Cellules (µL)</th><th class="l">Méthode</th></tr></thead><tbody>';
  let pbsa=0,tot=0;
  const sl=sampleList(), zv=+S.zeroVol;
  const nZeroLib=sl.filter(x=>x.strain==='banque'&&x.prim==='—').length;
  const nZeroCtl=sl.filter(x=>x.strain!=='banque'&&x.prim==='—').length;
  if(nZeroLib){ const v=nZeroLib*zv; tot+=v; pbsa+=v;
    h+='<tr><td>—</td><td>0</td><td>'+v.toFixed(0)+' <span class="muted">('+nZeroLib+' × '+zv+')</span></td>'+
      '<td class="l">'+vessel(v)+'</td><td class="l muted">—</td><td class="muted">—</td>'+
      '<td>'+v.toFixed(0)+'</td><td class="muted">—</td>'+
      '<td class="l muted">culot repris dans le PBSA</td></tr>'; }
  ladder().forEach(l=>{
    const r=rowFor(l), p=primaryTotal(r), src=aceSource(p.total,r.A);
    const m=addMode(p.total, r.cv);
    const rest=p.total-src.addUL-(m.direct?r.cv:0);
    pbsa+=rest; tot+=p.total;
    h+='<tr><td>'+l.toFixed(2)+'</td><td>'+sci(r.A)+'</td><td>'+p.total.toFixed(0)+'</td>'+
      '<td class="l">'+vessel(p.total)+'</td>'+
      '<td class="l">'+src.label+'</td>'+
      '<td'+(src.addUL<(+S.minPip)?' class="badcell bad"':'')+'>'+
        src.addUL.toFixed(src.addUL<10?2:1)+'</td>'+
      '<td>'+rest.toFixed(rest<100?1:0)+'</td>'+
      '<td>'+r.cv+'</td><td class="l'+(m.direct?'':' muted')+'">'+m.label+'</td></tr>';
  });
  const rx=controlRxns();
  if(rx.length){
    h+='<tr class="sub"><td class="l" colspan="9">Contrôles clonaux — solutions préparées à part</td></tr>';
    if(nZeroCtl){ const v=nZeroCtl*zv; tot+=v; pbsa+=v;
      h+='<tr><td>—</td><td>0</td><td>'+v.toFixed(0)+' <span class="muted">('+nZeroCtl+' × '+zv+')</span></td>'+
        '<td class="l">'+vessel(v)+'</td><td class="l muted">—</td><td class="muted">—</td>'+
        '<td>'+v.toFixed(0)+'</td><td class="muted">—</td>'+
        '<td class="l muted">culot repris dans le PBSA</td></tr>'; }
    rx.forEach(c=>{
      const cr=ctrlRow(c), src=aceSource(cr.total,cr.A);
      const m=addMode(c.unit, c.cv);
      const rest=cr.total-src.addUL-(m.direct?c.cv*c.n:0);
      tot+=cr.total; pbsa+=rest;
      h+='<tr><td>'+c.l.toFixed(2)+'</td><td>'+sci(cr.A)+'</td>'+
        '<td>'+cr.total.toFixed(0)+(c.n>1?' <span class="muted">('+c.n+' × '+c.unit+')</span>':'')+'</td>'+
        '<td class="l">'+vessel(cr.total)+'</td><td class="l">'+src.label+'</td>'+
        '<td'+(src.addUL<(+S.minPip)?' class="badcell bad"':'')+'>'+
          src.addUL.toFixed(src.addUL<10?2:1)+'</td>'+
        '<td>'+rest.toFixed(rest<100?1:0)+'</td>'+
        '<td>'+c.cv+(c.n>1?' <span class="muted">× '+c.n+'</span>':'')+'</td>'+
        '<td class="l'+(m.direct?'':' muted')+'">'+m.label+'</td></tr>';
    });
  }
  h+='</tbody><tfoot><tr><td class="l" colspan="2">Totaux</td><td>'+vol(tot)+'</td>'+
     '<td colspan="3"></td><td>'+vol(pbsa)+'</td><td colspan="2"></td></tr></tfoot>';
  t.innerHTML=h;

  const st=stockSummary();
  const neat=st.neat, d10=st.d10, d100=st.d100;
  const have=+S.stockHave;
  const mk=x=>x>0?Math.max(20,Math.ceil(x*2/10)*10):0;
  const make100=mk(d100);
  const make10=mk(d10+make100/10);
  const grand=neat+make10/10;
  const tight=grand>have*0.85;
  const short=grand>have;
  $('primeNotes').innerHTML=
    '<div class="twocol">'+
    '<div class="panel"><h3>Dilutions en série à préparer d’abord</h3>'+
	(d10>0?kv(cap(LIG())+' 1:10, préparer une quantité de ', vol(make10)+' ('+(make10/10).toFixed(1)+' µL de stock + '+
      (make10*0.9).toFixed(0)+' µL de PBSA)'):'')+
	(d100>0?kv(cap(LIG())+' 1:100, préparer une quantité de ', vol(make100)+' ('+(make100/10).toFixed(1)+' µL du 1:10 + '+
      (make100*0.9).toFixed(0)+' µL de PBSA)'):'')+
    kv('Utilisé directement dans les réactions', neat.toFixed(2)+' µL équivalent pur')+
    kv('Stock consommé au total','<span class="'+(short?'bad':(tight?'warn':'good'))+'">'+grand.toFixed(1)+
      ' µL = '+mass(massUG(grand))+'</span>')+
    kv('Vous avez déclaré disposer de', have.toFixed(0)+' µL = '+mass(massUG(have)))+
    '</div>'
    // '<div class="panel"><h3>Mode opératoire</h3>'+
    // '<p style="font-size:13.5px;color:var(--ink-2);margin:0">PBSA-BL mesuré dans le récipient, '+LIG()+
    // ' pipeté dedans, homogénéisation par retournements — vingt fois pour les coniques de 500 mL. '+
    // 'Dilutions dosées au nanodrop si possible.</p></div></div>';
}

function renderControls(){
  const rx=controlRxns(), cl=clones();
  const tag=$('nCtrl');
  if(tag) tag.textContent = rx.length+' réaction'+(rx.length>1?'s':'')+', '+cl.length+' souche'+(cl.length>1?'s':'');
  const picked=S.cloneConc;
  $('cloneConcPicker').innerHTML = ladder().map(l=>{
    const k=key(l), on=picked.includes(k);
    return '<button type="button" class="chip" data-cc="'+k+'" style="cursor:pointer;background:'+
      (on?'var(--teal-soft)':'transparent')+';'+(on?'border-color:var(--teal);color:var(--teal)':'')+
      '">e'+l+'</button>';
  }).join('');
  const nConc=S.cloneConc.filter(k=>ladder().map(key).includes(k)).length;
  const cellsNeeded=cl.length?(1+nConc)*(+S.cellVolCtl):0;
  $('ctrlNotes').innerHTML='';
}
function renderCellPrep(){
  const sl=sampleList(), cl=clones();
  const cultures=[{id:'__lib', nom:(S.libName||'LIB')+' (banque)',
                   need:sl.filter(x=>x.strain==='banque').reduce((a,b)=>a+b.cv,0),
                   resus:+S.libResus, tubes:2}]
    .concat(cl.map(c=>({id:c, nom:c,
                        need:sl.filter(x=>x.strain===c).reduce((a,b)=>a+b.cv,0),
                        resus:+S.cloneResus, tubes:1})));
  const dispo=c=>c.tubes*(+S.indVol)*1000;
  let h='', alerts=[];
  cultures.forEach(c=>{
    const od=S.odPrep[c.id]!=null?+S.odPrep[c.id]:'';
    const real=od?od*10:null;
    const pellet=real?c.resus/real:null;
    const spin=pellet&&pellet>1400?'15 mL, 3000 g, 5 min':'1,5 mL, 14000 g, 1 min';
    const manque=pellet!=null&&pellet>dispo(c);
    h+='<div class="cellrow"><div class="nm">'+c.nom+'</div>'+
      '<div class="od">DO 1:10 <input type="number" step="0.01" data-od="'+c.id+'" value="'+od+'"></div>'+
      (pellet!=null
        ? '<div class="out">culotter <b'+(manque?' class="bad"':'')+'>'+pellet.toFixed(0)+' µL</b>'+
          ' · '+spin+' · reprendre dans <b>'+vol(c.resus)+'</b></div>'
        : '<div class="out wait">reprise dans '+vol(c.resus)+', volume à culotter en attente de la DO</div>')+
      '</div>';
    if(manque) alerts.push('<span class="bad">'+c.nom+' : il faudrait culotter '+vol(pellet)+
      ' alors que la culture n\'en contient que '+vol(dispo(c))+'. La DO est trop basse pour atteindre '+
      vol(c.resus)+' à DO 1.</span>');
    if(c.resus<c.need) alerts.push('<span class="bad">'+c.nom+' : les échantillons demandent '+
      vol(c.need)+' de cellules, la reprise est fixée à '+vol(c.resus)+'.</span>');
    else if(c.resus<c.need*1.2) alerts.push('<span class="warn">'+c.nom+' : '+vol(c.need)+' nécessaires '+
      'pour '+vol(c.resus)+' de reprise, les lavages laisseront peu de marge.</span>');
  });
  $('cellPrep').innerHTML=h;
  $('cellPrepNote').innerHTML=alerts.length
    ? '<div class="note">'+alerts.join('<br>')+'</div>' : '';
  $('resusLine').textContent='('+vol(+S.libResus)+' pour la banque, '+vol(+S.cloneResus)+' par clone)';
}

function renderSecondary(){
  const sl=sampleList();
  const groups=[
    {n:'strep-PE seul', match:s=>s.sec==='strep-PE', pe:true,myc:false,fitc:false},
    {n:'myc-FITC seul', match:s=>s.sec==='myc-FITC', pe:false,myc:true,fitc:false},
    {n:'strep-PE + myc-FITC', match:s=>s.sec==='strep-PE + myc-FITC', pe:true,myc:true,fitc:false},
    {n:'FITC seul', match:s=>s.sec==='FITC', pe:false,myc:false,fitc:true},
    {n:'PBSA-BL (sans secondaire)', match:s=>s.sec==='—', pe:false,myc:false,fitc:false}
  ];
  let h='<thead><tr><th class="l">Secondaire</th><th>Échantillons</th><th>Volume + excédent (µL)</th>'+
    '<th>strep-PE (µL)</th><th>myc-FITC (µL)</th><th>FITC (µL)</th><th>PBSA-BL (µL)</th>'+
    '<th class="l">Récipient</th></tr></thead><tbody>';
  let tPE=0,tMyc=0,tFITC=0;
  groups.forEach(g=>{
    const n=sl.filter(g.match).length; if(!n)return;
    const v=n*(+S.secVol)*(1+(+S.secOver)/100);
    const pe=g.pe?v/(+S.dilPE):0, myc=g.myc?v/(+S.dilMyc):0, fi=g.fitc?v/(+S.dilFITC):0;
    tPE+=pe;tMyc+=myc;tFITC+=fi;
    h+='<tr><td class="l">'+g.n+'</td><td>'+n+'</td><td>'+v.toFixed(0)+'</td>'+
      '<td>'+(pe?pe.toFixed(2):'<span class="muted">—</span>')+'</td>'+
      '<td>'+(myc?myc.toFixed(2):'<span class="muted">—</span>')+'</td>'+
      '<td>'+(fi?fi.toFixed(2):'<span class="muted">—</span>')+'</td>'+
      '<td>'+(v-pe-myc-fi).toFixed(0)+'</td><td class="l">'+vessel(v)+'</td></tr>';
  });
  h+='</tbody><tfoot><tr><td class="l" colspan="3">Anticorps consommé</td><td>'+tPE.toFixed(1)+'</td>'+
     '<td>'+tMyc.toFixed(1)+'</td><td>'+tFITC.toFixed(1)+'</td><td></td><td></td></tr></tfoot>';
  $('secTable').innerHTML=h;
  window.__secUse={pe:tPE,myc:tMyc,fitc:tFITC};
}

function sampleRows(groups){
  let h='<thead><tr><th>#</th><th class="l">Échantillon</th><th class="l">Souche</th><th class="l">Primaire</th>'+
    '<th class="l">Secondaire</th><th>Cellules (µL)</th><th>Vol. marquage (µL)</th><th class="l">Récipient</th>'+
    '<th class="l">Rôle</th></tr></thead><tbody>';
  let i=0;
  groups.forEach(g=>{
    h+='<tr class="sub"><td class="l" colspan="9">'+g.title+'</td></tr>';
    g.rows.forEach(s=>{ i++;
    h+='<tr><td>'+i+'</td><td class="l">'+s.n+'</td><td class="l">'+s.strain+'</td>'+
      '<td class="l">'+s.prim+'</td><td class="l">'+s.sec+'</td><td>'+s.cv+'</td>'+
      '<td>'+s.tv.toFixed(0)+'</td><td class="l">'+vessel(s.tv)+'</td>'+
      '<td class="l muted">'+s.why+'</td></tr>';
    });
  });
  return h+'</tbody>';
}
function renderSamples(){
  const sl=orderedSamples();
  const sle=$('seedLabel');
  if(sle) sle.textContent='tirage n° '+(+S.seed||1);
  const flow=sl.filter(s=>s.strain==='banque' && s.conc==null);
  const lib=sl.filter(s=>s.strain==='banque' && s.conc!=null);
  const ctrl=sl.filter(s=>s.strain!=='banque');
  $('nSamples').textContent=sl.length+' échantillons';
  const libTot=sl.filter(s=>s.strain==='banque').reduce((a,b)=>a+b.cv,0);
  $('sampleTable').innerHTML=sampleRows([
      {title:'Contrôles de cytométrie — à passer en premier', rows:flow},
      {title:'Contrôles clonaux', rows:ctrl},
    {title:'Banque — ordre tiré au sort', rows:lib}

  ])+'<tfoot><tr><td colspan="5" class="l">Cellules de banque requises</td>'+
    '<td>'+libTot+'</td><td colspan="3"></td></tr></tfoot>';
}

/* volume of every primary solution to make, library and controls */
function allVols(){
  return ladder().map(l=>primaryTotal(rowFor(l)).total)
    .concat(controlRxns().map(c=>ctrlRow(c).total));
}
function renderPrep(){
  const sl=sampleList();
  const nRec=ladder().length*(+S.nPEbins)+(+S.nFITCbins)+1;   // recovery cultures
  const sec=window.__secUse||{pe:0,myc:0,fitc:0};
  const primTot=ladder().reduce((a,l)=>a+primaryTotal(rowFor(l)).total,0);
  const washes=sl.length*4000 + sl.length*500 + 20000;
  const pbsaTot=(primTot+washes)*1.6;
  const collect=(nRec+8)*5;                    // mL de 2x SDCAA + 1% BSA
  const ypd=collect*0.749, bsa=collect*0.25;
  const sdcaa=(nRec*(+S.outVol)) + 5*(1+clones().length) + 50;
  const sgdcaa=(+S.indVol)*(2+clones().length) + 4.5*(1+clones().length) + 10;
  const platesPre=String(S.plateSamples).split(',').filter(Boolean).length;
  const platesPost=Object.keys(S.plateFlags).filter(k=>S.plateFlags[k]).length||6;
  const nPlates=(platesPre+platesPost)*2;

  $('prepWare').innerHTML=[
    ['',nRec+8+' tubes à culture'],
      ['', (nRec+16)+' tubes FACS 5 mL polypropylène'],
      ['', (sl.length*3+20)+' tubes de 1,5 mL'],
      ['', allVols().filter(v=>v>50000).length+' × coniques 500 mL'],
      ['', allVols().filter(v=>v>1500&&v<=15000).length+' tubes de 15 mL'],
      ['', allVols().filter(v=>v>15000&&v<=50000).length+' tubes de 50 mL'],
      ['', (sl.length+10)+' cryotubes'],
  ].map(checkrow).join('');

  $('prepReag').innerHTML=[
    [vol(sec.myc*1.3),'myc-FITC'],
    [vol(sec.pe*1.3),'strep-RPE'],
    [vol(sec.fitc*1.3),'FITC libre'],
    [vol(+S.stockHave),cap(LIG())],
    [70 * sl.length/1000 + ' mL','glycérol 40 % pour les cryotubes']
  ].map(checkrow).join('');

  $('prepMedia').innerHTML=(S.doPlating?[
    [nPlates+' + 30 de réserve','boîtes SDCAA + Amp, coulées'],
    [nPlates+' + 30 de réserve','boîtes YPD + Amp, coulées']]:[]).concat([
    [(pbsaTot/1e6).toFixed(1)+' L','PBSA-BL'],
    [Math.ceil(sdcaa)+' mL','SDCAA-BL'],
    [Math.ceil(sgdcaa)+' mL','SGDCAA-BL'],
    [Math.ceil(ypd)+' mL','SDCAA 2,67×'],
    [Math.ceil(bsa)+' mL','BSA 4 % — 1 % final dans les tubes de collecte'],
    [Math.ceil(collect)+' mL','2× SDCAA + 1 % BSA, mélangés, gardés au froid'],
  ]).map(checkrow).join('');


  $('tubeCount').innerHTML=
    '<div class="grid">'+
    kv('Tubes de cellules (1,5 mL)',sl.length)+
    kv('Récipients de solution primaire',ladder().length+' banque + '+controlRxns().length+' contrôles')+
    kv('Tubes de collecte du tri',(ladder().length*(+S.nPEbins))+(+S.nFITCbins))+
    kv('Tubes à culture pour la récupération',nRec)+
    kv('Tubes de congélation',sl.length+' cryo + '+sl.length+' × 1,5 mL')+
    '</div>';
}
function checkrow(a){ return '<label class="checkrow"><input type="checkbox">'+
  '<span class="nm">'+a[1]+'</span><span class="qty">'+a[0]+'</span></label>'; }

function renderBackdil(){
  const sat=+S.satOD, tgt=+S.indOD, v=+S.indVol;
  const take=sat>0?(v*1000*tgt/sat):null;
  $('backdil').innerHTML='<div class="grid">'+
    kv('Prélever dans '+v+' mL de SGDCAA-BL', take?vol(take):'—')+
    kv('Cellules à l’induction', cells(v*(+S.cellsPerOD)*tgt)+' par tube')+
    kv('Au-dessus de la diversité', (v*(+S.cellsPerOD)*tgt/(+S.diversity)).toFixed(0)+'×')+
    '</div>';
}

function renderPrePlate(){
  const list=String(S.plateSamples).split(',').map(s=>s.trim()).filter(Boolean);
  let h='<thead><tr><th class="l">Échantillon</th><th>Cellules dans la réaction</th><th>Vol. marquage (µL)</th>'+
    '<th>Cellules par µL</th><th>Étaler (µL)</th><th class="l">Remarque</th></tr></thead><tbody>';
  list.forEach(k=>{
    const l=parseFloat(k); if(isNaN(l))return;
    const r=rowFor(l);
    const perUL=r.nCells/r.tv;
    const need=(+S.targetCol)/perUL;
    let note='directement depuis le tube';
    let shown=need;
    if(need<2){ shown=need*50; note='après une dilution au 1:50 — 2 µL dans 98 µL d’eau'; }
    else if(need>150){ shown=150; note='plafonné à 150 µL ; attendez-vous à '+(150*perUL).toFixed(0)+' colonies'; }
    h+='<tr><td class="l">'+(S.libName||'LIB')+'_'+l+'</td><td>'+cells(r.nCells)+'</td><td>'+r.tv.toFixed(0)+'</td>'+
      '<td>'+perUL.toFixed(perUL<10?2:0)+'</td><td>'+shown.toFixed(shown<10?1:0)+'</td>'+
      '<td class="l muted">'+note+'</td></tr>';
  });
  h+='</tbody>';
  $('preplate').innerHTML=h;
}

function sortRows(){
  const rows=[];
  rows.push({id:'FITC', label:(S.libName||'LIB')+'_mycFITC — par passage (pairs, puis impairs)',
    bins:Math.max(1,Math.round((+S.nFITCbins)/2))});
  runOrder().forEach(l=>rows.push({id:'e'+l,label:(S.libName||'LIB')+'_'+l,bins:+S.nPEbins}));
  return rows;
}
function renderSortLog(){
  const b=budget();
  $('fitcTarget').textContent=cells(b.sorted*1.6);
  $('peTarget').textContent=cells(+S.sortTarget);
  let h='<thead><tr><th class="l">Échantillon</th><th>Heure</th><th>Objectif</th><th>Réellement trié</th>'+
    '<th>Par bin</th><th>Récupéré par bin</th><th>× la diversité</th></tr></thead><tbody>';
  sortRows().forEach(r=>{
    const tgt=r.id==='FITC'?b.sorted*1.6:(+S.sortTarget);
    const act=S.sorted[r.id]!=null?+S.sorted[r.id]:'';
    const perBin=act?act/r.bins:null;
    const rec=perBin?perBin*(1-(+S.lossRecov)/100):null;
    const cov=rec?rec/(+S.diversity):null;
    const cls=cov==null?'':(cov<5?'badcell bad':(cov<10?'warncell warn':''));
    const tm=S.sortTime[r.id]!=null?S.sortTime[r.id]:'';
    h+='<tr><td class="l">'+r.label+'</td>'+
      '<td><input class="mini" type="time" data-tm="'+r.id+'" value="'+tm+'"></td>'+
      '<td class="muted">'+cells(tgt)+'</td>'+
      '<td><input class="mini wide" type="number" step="10000" data-sort="'+r.id+'" value="'+act+'"></td>'+
      '<td>'+(perBin?cells(perBin):'—')+'</td><td>'+(rec?cells(rec):'—')+'</td>'+
      '<td class="'+cls+'">'+(cov?cov.toFixed(1)+'×':'—')+'</td></tr>';
  });
  h+='</tbody>';
  $('sortLog').innerHTML=h;
}

function renderPostPlate(){
  const out=(+S.outVol)*1000;
  const qc=S.qcConc;
  const rows=[];
  for(let i=1;i<=(+S.nFITCbins);i++) rows.push({id:'FITC'+i, src:'FITC', bin:i, bins:+S.nFITCbins});
  for(let i=1;i<=(+S.nPEbins);i++) rows.push({id:qc+'_PE'+i, src:'e'+qc, bin:i, bins:+S.nPEbins});
  let h='<thead><tr><th class="l">Bin</th><th>Étaler ?</th><th>Cellules dans le bin</th><th>Par µL après croissance</th>'+
    '<th>Étaler (µL)</th><th class="l">Eau d’étalement</th></tr></thead><tbody>';
  rows.forEach(r=>{
    const tot=S.sorted[r.src]!=null?+S.sorted[r.src]/r.bins:null;
    const cellsIn=S.binCells[r.id]!=null?+S.binCells[r.id]:(tot||'');
    const perUL=cellsIn?cellsIn/out:null;
    const need=perUL?(+S.targetCol2)/perUL:null;
    const flag=S.plateFlags[r.id];
    h+='<tr><td class="l">'+r.id+'</td>'+
      '<td><input type="checkbox" data-pf="'+r.id+'"'+(flag?' checked':'')+'></td>'+
      '<td><input class="mini wide" type="number" step="10000" data-bc="'+r.id+'" value="'+cellsIn+'"></td>'+
      '<td>'+(perUL?perUL.toFixed(0):'—')+'</td>'+
      '<td>'+(need?need.toFixed(0):'—')+'</td>'+
      '<td class="l muted">'+(need?(need<50?'compléter à 50–150 µL avec de l’eau':'étaler tel quel'):'—')+'</td></tr>';
  });
  h+='</tbody>';
  $('postplate').innerHTML=h;
  const nRec=ladder().length*(+S.nPEbins)+(+S.nFITCbins)+1;
  $('nCulture').textContent=nRec;
  $('halfOut').textContent=vol(out/2);
  const uns=(400000/(+S.cellsPerOD))*1000;
  $('unsortedVol').textContent=vol(uns);
  $('nFreeze').textContent=nRec;
}

function renderODLog(){
  const LB=S.libName||'LIB';
  const rows=ladder().map(l=>LB+'_'+l).concat([LB+'_0']);
  for(let i=0;i<Math.ceil((+S.nFITCbins)/(+S.nPEbins));i++) rows.push('groupe FITC '+(i+1));
  rows.push('non trié');
  let h='<thead><tr><th class="l">Culture</th>';
  for(let i=1;i<=(+S.nPEbins);i++) h+='<th>bin '+i+'</th>';
  h+='</tr></thead><tbody>';
  rows.forEach(r=>{
    h+='<tr><td class="l">'+r+'</td>';
    for(let i=1;i<=(+S.nPEbins);i++){
      const k=r+'|'+i; const v=S.odLog[k]!=null?S.odLog[k]:'';
      const bad=v!==''&&(+v<0.9||+v>2);
      h+='<td class="'+(bad?'warncell':'')+'"><input class="mini" type="number" step="0.01" data-odl="'+k+'" value="'+v+'"></td>';
    }
    h+='</tr>';
  });
  h+='</tbody>';
  $('odLog').innerHTML=h;
}

function renderCounts(){
  let h='<thead><tr><th class="l">Échantillon</th><th>YPD-Amp</th><th>SDCAA</th><th>Attendu</th>'+
    '<th>% de l’attendu</th><th>% de perte de plasmide</th><th></th></tr></thead><tbody>';
  S.counts.forEach((c,i)=>{
    const y=+c.y, d=+c.d, e=+c.e;
    const pct=(d&&e)?100*d/e:null;
    const loss=(y&&d)?100*(1-d/y):null;
    const lcls=loss!=null&&loss>40?'warn':'';
    h+='<tr><td class="l"><input type="text" data-cs="'+i+'" value="'+(c.s||'')+'" style="width:190px;margin:0;padding:3px 6px;font-size:13px;font-family:var(--sans)"></td>'+
      '<td><input class="mini" type="number" data-cy="'+i+'" value="'+(c.y||'')+'"></td>'+
      '<td><input class="mini" type="number" data-cd="'+i+'" value="'+(c.d||'')+'"></td>'+
      '<td><input class="mini" type="number" data-ce="'+i+'" value="'+(c.e||'')+'"></td>'+
      '<td>'+(pct!=null?pct.toFixed(1)+'%':'—')+'</td>'+
      '<td class="'+lcls+'">'+(loss!=null?loss.toFixed(1)+'%':'—')+'</td>'+
      '<td><button class="iconbtn" data-cdel="'+i+'">×</button></td></tr>';
  });
  h+='</tbody>';
  $('countTable').innerHTML=h;
}

/* ---------------- wiring ---------------- */
const FIELDS=['anchor','indH','primH','statH','ligName','dispName','libName','mw','stockMode','stockVal',
  'stockHave','cellsPerOD','recPerCell','worstLogKd',
  'maxDepl','minPip','zeroVol','minHeadroom','diversity','cellVol','cellVolCtl','sortTarget','lossLabel','lossWaste','lossRecov',
  'nPEbins','logStart','logEnd','logStep','maxVessel','cloneNames','satOD','indOD','indVol',
  'secVol','secOver','dilPE','dilMyc','dilFITC','targetCol','plateSamples','outVol','targetCol2'];


/* ---------------- partage par URL ----------------
   Seuls les paramètres du protocole voyagent dans le lien. Les journaux du jour
   (tri, DO, colonies, heures) restent dans le navigateur de chacun.
   Les valeurs égales aux valeurs de référence sont omises pour garder l'URL courte. */
const SHARE_SCALARS = FIELDS.concat(['doPlating','seed','libResus','cloneResus']);
const SHARE_MAPS = ['vols','cellVols','cloneVols'];

function encMap(o){ return Object.entries(o||{}).map(([k,v])=>k+':'+v).join(','); }
function decMap(str){
  const m={};
  str.split(',').forEach(pair=>{
    const i=pair.lastIndexOf(':'); if(i<1) return;
    const n=Number(pair.slice(i+1)); if(isFinite(n)) m[pair.slice(0,i)]=n;
  });
  return m;
}

function shareURL(){
  const p=new URLSearchParams();
  SHARE_SCALARS.forEach(f=>{
    const v=S[f];
    if(v===''||v==null||String(v)===String(DEF[f])) return;
    p.set(f, typeof v==='boolean' ? (v?'1':'0') : String(v));
  });
  SHARE_MAPS.forEach(f=>{ if(encMap(S[f])!==encMap(DEF[f])) p.set(f, encMap(S[f])); });
  if(S.cloneConc.join(',')!==DEF.cloneConc.join(',')) p.set('cloneConc', S.cloneConc.join(','));
  p.set('plateFlags', Object.keys(S.plateFlags).filter(k=>S.plateFlags[k]).join(','));
  return location.href.split(/[?#]/)[0]+'?'+p.toString();
}

function readShare(){
  const p=new URLSearchParams(location.search);
  if(![...p.keys()].some(k=>k in DEF)) return null;
  const o={};
  SHARE_SCALARS.forEach(f=>{
    if(!p.has(f)) return;
    const raw=p.get(f), d=DEF[f];
    if(typeof d==='boolean') o[f]= raw==='1'||raw==='true';
    else if(typeof d==='number'){ const n=Number(raw); if(raw!==''&&isFinite(n)) o[f]=n; }
    else o[f]=raw;
  });
  SHARE_MAPS.forEach(f=>{ if(p.has(f)) o[f]=decMap(p.get(f)); });
  if(p.has('cloneConc')) o.cloneConc=p.get('cloneConc').split(',').map(x=>x.trim()).filter(Boolean);
  if(p.has('plateFlags')){ o.plateFlags={};
    p.get('plateFlags').split(',').filter(Boolean).forEach(k=>o.plateFlags[k]=true); }
  return o;
}

function pushToInputs(){
  FIELDS.forEach(f=>{ const n=$(f); if(n && S[f]!=null) n.value=S[f]; });
  const dp=$('doPlating'); if(dp) dp.checked=!!S.doPlating;
  // $('indH2').textContent=S.indH; $('indH3').textContent=S.indH;
  // $('primH2').textContent=S.primH; $('statH2').textContent=S.statH;
  $('cloneList2').textContent='('+clones().join(', ')+')';
  $('c-cells').textContent=sci(+S.cellsPerOD,1);
  $('c-rec').textContent=sci(+S.recPerCell,1);
  $('thawVol').textContent='150 µL';
}

function captureFocus(){
  const a=document.activeElement;
  if(!a||a.tagName!=='INPUT')return null;
  let sel=null;
  if(a.id) sel='#'+a.id;
  else { const k=Object.keys(a.dataset)[0];
    if(k) sel='[data-'+k.replace(/[A-Z]/g,m=>'-'+m.toLowerCase())+'="'+a.dataset[k]+'"]'; }
  if(!sel)return null;
  let s=null,e=null; try{ s=a.selectionStart; e=a.selectionEnd; }catch(_){}
  return {sel,raw:a.value,s,e};
}
function restoreFocus(f){
  if(!f)return; let n=null; try{ n=document.querySelector(f.sel); }catch(_){}
  if(!n)return;
  if(n.value!==f.raw) n.value=f.raw;
  n.focus();
  try{ if(f.s!=null) n.setSelectionRange(f.s,f.e); }catch(_){}
}
function applyPlating(){
  const on=!!S.doPlating;
  document.querySelectorAll('.plating').forEach(n=>{ n.style.display = on?'':'none'; });
  const t=$('recovTitle');
  if(t) t.textContent = on ? 'Récupération et étalement' : 'Récupération';
  const m=$('mediaTitle');
  if(m) m.textContent = on ? 'Milieux et boîtes' : 'Milieux';
}
function applyNames(){
  document.querySelectorAll('.lig').forEach(n=>n.textContent=LIG());
  document.querySelectorAll('.disp').forEach(n=>n.textContent=DISP());
  const lbl=$('stockValLabel');
  if(lbl) lbl.textContent = S.stockMode==='mass'
    ? 'Concentration du stock (mg/mL)' : 'Concentration du stock (µM)';
  const dv=$('stockDerived');
  if(dv) dv.innerHTML='<div class="grid" style="margin-top:14px">'+
    kv('Molarité du stock', stockUM().toFixed(3)+' µM')+
    kv('Concentration massique', stockMgMl().toFixed(3)+' mg/mL')+
    kv('Dans le tube', (+S.stockHave).toFixed(0)+' µL = '+mass(massUG(+S.stockHave)))+
    '</div>';
}
function renderAll(){
  const f=captureFocus();
  applyNames(); applyPlating();
  renderClock(); renderBudget(); renderCurve(); renderPrimary(); renderCellPrep();
  renderControls(); renderSecondary(); renderSamples(); renderPrep(); renderBackdil(); renderPrePlate();
  renderSortLog(); renderPostPlate(); renderODLog(); renderCounts();
  restoreFocus(f);
  save();
}

document.addEventListener('input', e=>{
  const t=e.target, d=t.dataset;
  if(t.id && FIELDS.includes(t.id)){ S[t.id]= t.type==='number'? (t.value===''?'':+t.value) : t.value; renderAll(); return; }
  if(d.tv!=null){ S.vols[d.tv]=+t.value; renderAll(); }
  else if(d.cv!=null){ S.cellVols[d.cv]=+t.value; renderAll(); }
  else if(d.clv!=null){ S.cloneVols[d.clv]=+t.value; renderAll(); }
  else if(d.od!=null){ S.odPrep[d.od]=t.value===''?null:+t.value; renderAll(); }
  else if(d.res!=null){ if(d.res==='__lib') S.libResus=+t.value; else S.cloneResus=+t.value; renderAll(); }
  else if(d.sort!=null){ S.sorted[d.sort]=t.value===''?null:+t.value; renderAll(); }
  else if(d.tm!=null){ S.sortTime[d.tm]=t.value; save(); }
  else if(d.bc!=null){ S.binCells[d.bc]=t.value===''?null:+t.value; renderAll(); }
  else if(d.odl!=null){ S.odLog[d.odl]=t.value===''?'':+t.value; save(); }
  else if(d.cs!=null){ S.counts[+d.cs].s=t.value; save(); }
  else if(d.cy!=null){ S.counts[+d.cy].y=t.value; renderCounts(); save(); }
  else if(d.cd!=null){ S.counts[+d.cd].d=t.value; renderCounts(); save(); }
  else if(d.ce!=null){ S.counts[+d.ce].e=t.value; renderCounts(); save(); }
});
document.addEventListener('change', e=>{
  const d=e.target.dataset;

  if(d.pf!=null){ S.plateFlags[d.pf]=e.target.checked; save(); }
  if(e.target.id==='doPlating'){ S.doPlating=e.target.checked; renderAll(); }
});
document.addEventListener('click', e=>{
  const d=e.target.dataset;
  if(d.cc!=null){
    const k=d.cc;
    if(S.cloneConc.includes(k)) S.cloneConc=S.cloneConc.filter(x=>x!==k);
    else { S.cloneConc.push(k); if(!S.cloneVols[k]) S.cloneVols[k]=1000; }
    renderAll(); return;
  }
  if(d.cdel!=null){ S.counts.splice(+d.cdel,1); renderCounts(); save(); }
});

$('autoVol').addEventListener('click',()=>{
  ladder().forEach(l=>{ const k=key(l);
    const cv=S.cellVols[k]!=null?+S.cellVols[k]:+S.cellVol;
    S.vols[k]=autoVolume(l,cv); });
  renderAll();
});
$('restoreVol').addEventListener('click',()=>{
  S.vols=Object.assign({},DEFAULT_VOLS); S.cellVols=Object.assign({},DEFAULT_CELLS); renderAll();
});
$('reshuffle').addEventListener('click',()=>{
  S.seed=Math.floor(Math.random()*1e6)+1; renderAll();
});
$('addCount').addEventListener('click',()=>{ S.counts.push({s:'',y:'',d:'',e:''}); renderCounts(); save(); });
$('copySamples').addEventListener('click',()=>{
  const sl=orderedSamples();
  const tsv=['ordre\techantillon\tsouche\tprimaire\tsecondaire\tcellules_uL\tmarquage_uL\trecipient\trole']
    .concat(sl.map((s,i)=>[i+1,s.n,s.strain,s.prim,s.sec,s.cv,s.tv.toFixed(0),vessel(s.tv),s.why].join('\t'))).join('\n');
  navigator.clipboard.writeText(tsv).then(()=>{
    const b=$('copySamples'); const o=b.textContent; b.textContent='Copié';
    setTimeout(()=>b.textContent=o,1400); }).catch(()=>{});
});
$('shareBtn').addEventListener('click',()=>{
  const url=shareURL(), b=$('shareBtn'), o=b.textContent;
  const done=()=>{ b.textContent='Lien copié'; setTimeout(()=>b.textContent=o,1600); };
  if(navigator.clipboard) navigator.clipboard.writeText(url).then(done).catch(()=>prompt('Copier ce lien :',url));
  else prompt('Copier ce lien :',url);
});
$('printBtn').addEventListener('click',()=>window.print());
$('resetBtn').addEventListener('click',()=>{
  if(confirm('Effacer toutes les valeurs saisies sur cette page et revenir aux valeurs de référence ?')){
    S=JSON.parse(JSON.stringify(DEF)); try{localStorage.removeItem(STORE_KEY);}catch(e){}
    pushToInputs(); renderAll(); }
});
$('themeBtn').addEventListener('click',()=>{
  const r=document.documentElement;
  const dark=r.getAttribute('data-theme')==='dark' ||
    (!r.getAttribute('data-theme') && matchMedia('(prefers-color-scheme:dark)').matches);
  r.setAttribute('data-theme', dark?'light':'dark');
  $('themeBtn').textContent = dark?'Sombre':'Clair';
  try{ localStorage.setItem('titeseq-theme', dark?'light':'dark'); }catch(e){}
});

/* boot */
load();
{
  const shared=readShare();
  if(shared){
    let hadLocal=false; try{ hadLocal=!!localStorage.getItem(STORE_KEY); }catch(e){}
    if(!hadLocal || confirm('Ce lien contient un protocole. Remplacer les paramètres enregistrés dans ce navigateur ?\n\nLes journaux (tri, DO, colonies) seront remis à zéro.')){
      S=Object.assign(JSON.parse(JSON.stringify(DEF)), shared);
    }
    // on retire les paramètres de l'URL pour qu'un rechargement ne réécrase pas les modifications
    history.replaceState(null,'',location.pathname+location.hash);
  }
}
try{ const th=localStorage.getItem('titeseq-theme');
  if(th){ document.documentElement.setAttribute('data-theme',th);
    $('themeBtn').textContent = th==='dark'?'Clair':'Sombre'; } }catch(e){}
if(!S.anchor){
  const d=new Date(); d.setDate(d.getDate()+2); d.setHours(11,0,0,0);
  S.anchor=new Date(d.getTime()-d.getTimezoneOffset()*60000).toISOString().slice(0,16);
}
if(!Object.keys(S.plateFlags).length){
  ['FITC1','FITC3','FITC5','FITC7','9.5_PE1','9.5_PE4'].forEach(k=>S.plateFlags[k]=true);
}
pushToInputs();
renderAll();
