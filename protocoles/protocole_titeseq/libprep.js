/* ---------- helpers ---------- */
const $ = id => document.getElementById(id);
function vol(u){ if(u==null||!isFinite(u))return '—';
  if(u>=1000) return (u/1000).toFixed(u>=10000?0:2)+' mL';
  if(u>=10) return u.toFixed(0)+' µL';
  if(u>=1) return u.toFixed(2)+' µL';
  return u.toFixed(2)+' µL'; }
function cells(n){ if(!isFinite(n)||n==null)return '—';
  if(n>=1e9) return (n/1e9).toFixed(2)+'e9';
  if(n>=1e6) return (n/1e6).toFixed(n>=1e7?1:2)+'e6';
  if(n>=1e3) return (n/1e3).toFixed(0)+'e3';
  return n.toFixed(0); }
function kv(k,v){ return '<div class="kv"><span>'+k+'</span><span class="v">'+v+'</span></div>'; }
function checkrow(a){ return '<label class="checkrow"><input type="checkbox">'+
  '<span class="nm">'+a[1]+'</span><span class="qty">'+a[0]+'</span></label>'; }
function nf(x){ return (+x).toLocaleString('fr-FR'); }

/* ---------- state ---------- */
const DEF = {
  nRxn:null, over1:10, over2:18, dead:1,
  elut:18, tmpl:5, washVol:700, kitSize:50,
  p1Primer:1, p2Primer:2.5,
  cyc2:34, ngUl:20, plasmidBp:6000,
  waterRatio:1, bead1:1.2, elute1:35, recover1:33,
  gelLoad:5, bead2:0.85, elute2:50,
  totalReads:400, minPip:2, availUL:45, ampBp:400,
  poolVol:200, cutLow:0.55, cutHigh:0.85,
  finalNg:10, finalVol:30, corrFactor:1, phix:20,
  share:{}, conc:{}, binCells:{}
};
let S = JSON.parse(JSON.stringify(DEF));
let LINK = null;   /* state read from the sorting page */

function save(){ try{ localStorage.setItem('titeseq-libprep-fr-v1', JSON.stringify(S)); }catch(e){} }
function load(){ try{ const r=localStorage.getItem('titeseq-libprep-fr-v1');
  if(r) S=Object.assign(JSON.parse(JSON.stringify(DEF)), JSON.parse(r)); }catch(e){} }
function loadLink(){
  try{ const r=localStorage.getItem('titeseq-rbd-fr-v1'); if(!r) return null;
    const t=JSON.parse(r);
    const a=+t.logStart, b=+t.logEnd, st=Math.abs(+t.logStep)||1;
    const lad=[]; for(let x=a;x<=b+1e-9;x+=st) lad.push(Math.round(x*1000)/1000);
    return {lad, nPE:+t.nPEbins||4, nFITC:+t.nFITCbins||8, lib:t.libName||'LIB',
      sorted:t.sorted||{}, lossRecov:+t.lossRecov||30, diversity:+t.diversity||65536};
  }catch(e){ return null; }
}

/* ---------- the wells ---------- */
function bins(){
  const out=[];
  if(LINK){
    for(let i=1;i<=LINK.nFITC;i++) out.push({id:'FITC'+i, label:LINK.lib+'_FITC'+i, src:'FITC', bins:LINK.nFITC});
    LINK.lad.forEach(l=>{ for(let i=1;i<=LINK.nPE;i++)
      out.push({id:'e'+l+'_PE'+i, label:LINK.lib+'_'+l+'_PE'+i, src:'e'+l, bins:LINK.nPE}); });
    out.push({id:'unsorted', label:LINK.lib+'_non trié', src:null, bins:1});
  } else {
    const n=nRxn();
    for(let i=1;i<=n;i++) out.push({id:'w'+i, label:'puits '+i, src:null, bins:1});
  }
  return out;
}
function binCellCount(b){
  if(S.binCells[b.id]!=null) return +S.binCells[b.id];
  if(LINK && b.src && LINK.sorted[b.src]!=null)
    return (+LINK.sorted[b.src])/b.bins*(1-LINK.lossRecov/100);
  return null;
}
function nRxn(){
  if(S.nRxn) return +S.nRxn;
  if(LINK) return LINK.lad.length*LINK.nPE + LINK.nFITC + 1;
  return 33;
}
function nMM(){ return nRxn()*(1+(+S.over1)/100); }
function nMM2(){ return nRxn()*(1+(+S.over2)/100); }

/* ---------- render ---------- */
function renderLink(){
  $('linkState').innerHTML = LINK
    ? '<div class="callout">Échelle reprise de la page de tri : '+LINK.lad.length+
      ' concentrations × '+LINK.nPE+' bins PE, '+LINK.nFITC+' bins d\'expression, plus le non trié.</div>'
    : '<div class="callout w">La page de tri n\'a pas été trouvée dans ce navigateur. Saisissez le '+
      'nombre de réactions à la main ; tout le reste suivra.</div>';
}

function renderScale(){
  const n=nRxn(), strips=Math.ceil(n/8), plates=Math.ceil(n/96);
  $('nRxnTag').textContent=n+' réactions';
  $('scaleOut').innerHTML='<div class="grid" style="margin-top:14px">'+
    kv('Réactions', n)+
    kv('Barrettes de 8', strips)+
    kv('Plaques de 96', plates+(plates>1?' (attention au plan de plaque)':''))+
    kv('Équivalents MM, PCR 1', nMM().toFixed(1))+
    kv('Équivalents MM, PCR 2', nMM2().toFixed(1))+
    '</div>'+
    (n%8?'<div class="callout w">'+(8-n%8)+' puits vides dans la dernière barrette. Notez lesquels, '+
      'ils se remplissent tout seuls au moment du pooling.</div>':'');
}

function mmTable(el, rows, perRxn, factor, title){
  let h='<thead><tr><th class="l">Composant</th><th>Par réaction (µL)</th><th>× '+factor.toFixed(1)+
    ' (µL)</th></tr></thead><tbody>';
  let tot=0;
  rows.forEach(r=>{ tot+=r[1];
    h+='<tr><td class="l">'+r[0]+'</td><td>'+r[1].toFixed(1)+'</td><td>'+(r[1]*factor).toFixed(0)+
      '</td></tr>'; });
  h+='</tbody><tfoot><tr><td class="l">'+title+'</td><td>'+tot.toFixed(1)+'</td><td>'+
     (tot*factor).toFixed(0)+'</td></tr></tfoot>';
  $(el).innerHTML=h;
  return tot;
}

function renderPCR1(){
  const f=nMM();
  const mm=mmTable('mm1',[
    ['Tampon Q5 5×',4],['Eau',8.4],['dNTP',0.4],['Q5 polymérase',0.2]
  ],null,f,'Master mix par puits');
  $('mm1Txt').textContent=mm.toFixed(1)+' µL';
  $('tmplTxt').textContent=(+S.tmpl).toFixed(0)+' µL';
  return mm;
}
function renderPCR2(){
  const f=nMM2();
  const mm=mmTable('mm2',[
    ['Tampon KAPA 5×',10],['dNTP',1],['KAPA HiFi',1]
  ],null,f,'Master mix par puits');
  $('mm2Txt').textContent=mm.toFixed(1)+' µL';
  return mm;
}

function rxn1Vol(){ return 13 + (+S.tmpl) + 2*(+S.p1Primer); }
function rxn2Vol(){ return 12 + (+S.recover1) + 2*(+S.p2Primer); }

function renderPrimers(){
  const n=nRxn();
  let h='<thead><tr><th class="l">Étape</th><th>Par puits (µL)</th><th>Puits</th>'+
    '<th>Total par amorce colonne (µL)</th><th>Total par amorce ligne (µL)</th>'+
    '<th class="l">Remarque</th></tr></thead><tbody>';
  const perCol=Math.ceil(n/8), perRow=Math.min(n,8);
  [['PCR 1',+S.p1Primer],['PCR 2',+S.p2Primer]].forEach(([nm,v])=>{
    h+='<tr><td class="l">'+nm+'</td><td>'+v.toFixed(1)+'</td><td>'+n+'</td>'+
      '<td>'+(v*perRow).toFixed(1)+'</td><td>'+(v*perCol).toFixed(1)+'</td>'+
      '<td class="l muted">une amorce colonne et une amorce ligne par puits</td></tr>';
  });
  h+='</tbody>';
  $('primerTable').innerHTML=h;
}

function renderBead1(){
  const r=rxn1Vol(), water=r*(+S.waterRatio), after=r+water, beads=after*(+S.bead1);
  let h='<thead><tr><th class="l">Étape</th><th>Par puits (µL)</th><th>Total (µL)</th>'+
    '<th class="l">Remarque</th></tr></thead><tbody>';
  const n=nRxn();
  [['Volume de réaction après PCR 1',r,''],
   ['Eau ajoutée',water,'ratio 1:1'],
   ['Billes',beads,(+S.bead1).toFixed(2)+'× du volume dilué'],
   ['Lavage EtOH 85 %',150,'une fois, séchage jusqu\'à craquelure'],
   ['Élution',+S.elute1,'tampon Geneaid'],
   ['Récupéré',+S.recover1,'vers la plaque PCR 2']
  ].forEach(x=>{
    h+='<tr><td class="l">'+x[0]+'</td><td>'+x[1].toFixed(1)+'</td><td>'+((x[1]*n)/1000).toFixed(2)+
      ' mL</td><td class="l muted">'+x[2]+'</td></tr>';
  });
  h+='</tbody>';
  $('bead1Table').innerHTML=h;
}

function renderBead2(){
  const r=rxn2Vol(), afterGel=r-(+S.gelLoad);
  const n=nRxn();
  let h='<thead><tr><th class="l">Cas</th><th>Volume restant (µL)</th><th>Billes '+(+S.bead2).toFixed(2)+
    '× (µL)</th><th>Total billes (µL)</th></tr></thead><tbody>'+
    '<tr><td class="l">Gel après nettoyage</td><td>'+r.toFixed(1)+'</td><td>'+(r*(+S.bead2)).toFixed(2)+
      '</td><td>'+(r*(+S.bead2)*n).toFixed(0)+'</td></tr>'+
    '<tr><td class="l">Gel avant nettoyage ('+(+S.gelLoad)+' µL chargés)</td><td>'+afterGel.toFixed(1)+
      '</td><td>'+(afterGel*(+S.bead2)).toFixed(2)+'</td><td>'+(afterGel*(+S.bead2)*n).toFixed(0)+
      '</td></tr></tbody>';
  $('bead2Table').innerHTML=h;
  $('bead2Note').innerHTML='<div class="callout">Élution dans '+vol(+S.elute2)+' par puits, soit '+
    ((+S.elute2)*n/1000).toFixed(1)+' mL de tampon. Le volume de billes est recalculé à partir du volume '+
    'réellement restant : si vous chargez le gel avant, c\'est la deuxième ligne qui s\'applique.</div>';
}

function renderCycles(){
  const ng=+S.ngUl, ul=+S.tmpl, bp=+S.plasmidBp;
  const mols = ng*ul*1e-9/(660*bp) * 6.022e23;
  const dupl = LINK ? mols/(+LINK.diversity) : null;
  $('cycOut').innerHTML='<div class="grid" style="margin-top:14px">'+
    kv('Molécules entrant en PCR 1', cells(mols))+
    (dupl!=null?kv('Par variant de la banque','<span class="'+(dupl<10?'bad':(dupl<100?'warn':'good'))+'">'+
      dupl.toFixed(0)+'×</span>'):'')+
    kv('Cycles PCR 2', (+S.cyc2)+' (×'+Math.pow(2,+S.cyc2).toExponential(1)+' au maximum)')+
    '</div>'+
    (dupl!=null&&dupl<10
      ? '<div class="callout d">Moins de 10 copies par variant entrent en PCR. Le goulot est ici, pas '+
        'dans le séquençage : augmenter la profondeur ne récupérera rien. Augmentez le volume de matrice, '+
        'ou concentrez l\'éluat.</div>'
      : '<div class="callout">Un nombre de cycles fixe pousse les puits les plus concentrés en plateau, '+
        'où naissent les chimères entre variants. Si vous pouvez suivre quelques puits en qPCR, arrêtez '+
        'la PCR 2 en phase exponentielle plutôt qu\'à '+(+S.cyc2)+' cycles.</div>');
}

function renderPool(){
  const bs=bins();
  const counts=bs.map(binCellCount);
  const known=counts.filter(c=>c!=null&&c>0);
  const totalCells=known.reduce((a,b)=>a+b,0);
  /* read share: proportional to recovered cells unless overridden */
  const shares=bs.map((b,i)=>{
    if(S.share[b.id]!=null) return +S.share[b.id];
    if(counts[i]!=null&&totalCells>0) return 100*counts[i]/totalCells;
    return 100/bs.length;
  });
  const sumShare=shares.reduce((a,b)=>a+b,0)||1;
  /* volumes: mass in proportion to share, scaled so the smallest pipetted volume is respected */
  const concs=bs.map(b=>S.conc[b.id]!=null?+S.conc[b.id]:null);
  let K=0;
  bs.forEach((b,i)=>{ const c=concs[i]; if(c&&c>0){ const k=(+S.minPip)*c/(shares[i]/sumShare); if(k>K)K=k; } });
  let h='<thead><tr><th class="l">Fenêtre</th><th>Cellules récupérées</th><th>Part de lectures (%)</th>'+
    '<th>Lectures (M)</th><th>Lectures / cellule</th><th>Concentration (ng/µL)</th>'+
    '<th>À prélever (µL)</th></tr></thead><tbody>';
  let over=0, poolSum=0;
  bs.forEach((b,i)=>{
    const f=shares[i]/sumShare;
    const reads=f*(+S.totalReads);
    const perCell=counts[i]?reads*1e6/counts[i]:null;
    const c=concs[i];
    const v=(c&&c>0&&K>0)?K*f/c:null;
    if(v!=null){ poolSum+=v; if(v>(+S.availUL)) over++; }
    const vcls=v!=null&&v>(+S.availUL)?'badcell bad':'';
    const pcls=perCell!=null&&perCell<2?'warncell warn':'';
    h+='<tr><td class="l">'+b.label+'</td>'+
      '<td><input class="mini wide" type="number" step="10000" data-bc="'+b.id+'" value="'+
        (counts[i]!=null?Math.round(counts[i]):'')+'"></td>'+
      '<td><input class="mini" type="number" step="0.5" data-sh="'+b.id+'" value="'+
        (shares[i]).toFixed(2)+'"></td>'+
      '<td>'+reads.toFixed(1)+'</td>'+
      '<td class="'+pcls+'">'+(perCell!=null?perCell.toFixed(1):'—')+'</td>'+
      '<td><input class="mini" type="number" step="0.1" data-cc="'+b.id+'" value="'+
        (c!=null?c:'')+'"></td>'+
      '<td class="'+vcls+'">'+(v!=null?v.toFixed(1):'—')+'</td></tr>';
  });
  h+='</tbody><tfoot><tr><td class="l" colspan="6">Volume du pool</td><td>'+
     (poolSum?poolSum.toFixed(0)+' µL':'—')+'</td></tr></tfoot>';
  $('poolTable').innerHTML=h;
  $('poolNote').innerHTML=
    (over?'<div class="callout d">'+over+' fenêtre'+(over>1?'s demandent':' demande')+' plus que les '+
      vol(+S.availUL)+' disponibles. Baissez leur part, ou concentrez-les.</div>':'')+
    '<div class="callout">Les parts par défaut suivent les cellules récupérées. La colonne lectures par '+
    'cellule est le contrôle utile : en dessous de 2, la fenêtre est sous-échantillonnée ; très au-dessus '+
    'de 10, les lectures relisent les mêmes molécules et seraient mieux ailleurs.</div>';
}

function renderTwoSided(){
  const V=+S.poolVol, lo=+S.cutLow, hi=+S.cutHigh;
  const b1=V*lo, after=V+b1, b2=V*(hi-lo);
  $('twoSided').innerHTML='<div class="grid" style="margin-top:14px">'+
    kv('1 · billes à '+lo.toFixed(2)+'×', vol(b1)+' — garder le surnageant')+
    kv('Volume du surnageant', '≈ '+vol(after))+
    kv('2 · billes supplémentaires', vol(b2)+' — garder les billes')+
    kv('Billes au total', vol(b1+b2))+
    '</div>'+
    (hi<=lo?'<div class="callout d">La coupe haute doit être supérieure à la coupe basse.</div>':'');
}

function renderDilution(){
  const ng=+S.finalNg, bp=+S.ampBp, cf=+S.corrFactor||1;
  const nM = bp>0 ? (ng*1e6)/(660*bp) * cf : null;
  const target=10, Vf=+S.finalVol;
  const take = (nM&&nM>target) ? target*Vf/nM : null;
  $('dilOut').innerHTML='<div class="grid" style="margin-top:14px">'+
    kv('Molarité du pool', nM!=null?nM.toFixed(1)+' nM':'—')+
    kv('Pour '+Vf+' µL à 10 nM', take!=null?vol(take)+' de pool + '+vol(Vf-take)+' de tampon':'—')+
    kv('PhiX demandé', (+S.phix)+' %')+
    '</div>'+
    (nM!=null&&nM<target?'<div class="callout w">Le pool est déjà sous 10 nM. Soumettez-le tel quel en '+
      'indiquant sa molarité réelle.</div>':'');
}

function renderPrep(){
  const n=nRxn();
  const beads1=rxn1Vol()*(1+(+S.waterRatio))*(+S.bead1)*n;
  const beads2=rxn2Vol()*(+S.bead2)*n;
  const kits=Math.ceil(n/(+S.kitSize||50));
  $('prepWare').innerHTML=[
    [Math.ceil(n/8)+' barrettes','tubes PCR profil haut, bouchons plats'],
    [kits+' kit'+(kits>1?'s':''),'Zymo yeast miniprep II ('+(+S.kitSize)+' colonnes)'],
    ['1','plaque magnétique et multicanal'],
    [Math.ceil(n/96)+' plaque'+(Math.ceil(n/96)>1?'s':''),'plaques PCR pour les éluats'],
    ['1','gel 1 % au BET, coulé pendant la PCR 2']
  ].map(checkrow).join('');
  $('prepReag').innerHTML=[
    [vol(beads1+beads2),'billes A-line (plateforme), + marge'],
    [vol(((+S.elute1)+(+S.elute2))*n+(+S.elut)*n),'tampon d\'élution Geneaid'],
    [vol(150*n*2),'EtOH 85 % ou tampon W2'],
    [vol(13*nMM()),'master mix Q5 (voir PCR 1)'],
    [vol(12*nMM2()),'master mix KAPA (voir PCR 2)'],
    ['2 × '+n+' µL','amorces index, arrayées la veille']
  ].map(checkrow).join('');
}

function renderAll(){
  const a=document.activeElement;
  const sel=a&&a.tagName==='INPUT'?(a.id?'#'+a.id:(Object.keys(a.dataset)[0]?'[data-'+Object.keys(a.dataset)[0]+'="'+a.dataset[Object.keys(a.dataset)[0]]+'"]':null)):null;
  const raw=a&&a.tagName==='INPUT'?a.value:null;
  renderLink(); renderScale(); renderPrep(); renderPrimers();
  renderPCR1(); renderBead1(); renderPCR2(); renderCycles();
  renderBead2(); renderPool(); renderTwoSided(); renderDilution();
  $('washVolTxt').textContent=(+S.washVol)+' µL';
  $('elutTxt').textContent=(+S.elut)+' µL';
  $('cyc2Txt').textContent=(+S.cyc2);
  $('minioOut').innerHTML='<div class="grid" style="margin-top:14px">'+
    kv('Colonnes nécessaires', nRxn())+
    kv('Éluat par puits', vol(+S.elut))+
    kv('Engagé en PCR 1', vol(+S.tmpl)+' ('+(100*(+S.tmpl)/(+S.elut)).toFixed(0)+' % de l\'éluat)')+
    kv('Archivé à −20 °C', vol((+S.elut)-(+S.tmpl)))+
    '</div>';
  if(sel){ let n=null; try{ n=document.querySelector(sel); }catch(e){}
    if(n){ if(raw!=null&&n.value!==raw) n.value=raw; n.focus(); } }
  save();
}

/* ---------- wiring ---------- */
const FIELDS=['nRxn','over1','over2','dead','elut','tmpl','washVol','kitSize','cyc2','ngUl','plasmidBp',
  'totalReads','minPip','availUL','ampBp','poolVol','cutLow','cutHigh','finalNg','finalVol','corrFactor','phix'];
function pushToInputs(){ FIELDS.forEach(f=>{ const n=$(f); if(n&&S[f]!=null) n.value=S[f]; }); }

document.addEventListener('input', e=>{
  const t=e.target, d=t.dataset;
  if(t.id&&FIELDS.includes(t.id)){ S[t.id]= t.value===''?(t.id==='nRxn'?null:''):+t.value; renderAll(); return; }
  if(d.sh!=null){ S.share[d.sh]=t.value===''?null:+t.value; renderAll(); }
  else if(d.cc!=null){ S.conc[d.cc]=t.value===''?null:+t.value; renderAll(); }
  else if(d.bc!=null){ S.binCells[d.bc]=t.value===''?null:+t.value; renderAll(); }
});
$('printBtn').addEventListener('click',()=>window.print());
$('resetBtn').addEventListener('click',()=>{
  if(confirm('Effacer toutes les valeurs saisies sur cette page ?')){
    S=JSON.parse(JSON.stringify(DEF));
    try{ localStorage.removeItem('titeseq-libprep-fr-v1'); }catch(e){}
    pushToInputs(); renderAll(); }
});
$('themeBtn').addEventListener('click',()=>{
  const r=document.documentElement;
  const dark=r.getAttribute('data-theme')==='dark' ||
    (!r.getAttribute('data-theme')&&matchMedia('(prefers-color-scheme:dark)').matches);
  r.setAttribute('data-theme', dark?'light':'dark');
  $('themeBtn').textContent = dark?'Sombre':'Clair';
  try{ localStorage.setItem('titeseq-theme', dark?'light':'dark'); }catch(e){}
});

/* boot */
load();
LINK=loadLink();
try{ const th=localStorage.getItem('titeseq-theme');
  if(th){ document.documentElement.setAttribute('data-theme',th);
    $('themeBtn').textContent = th==='dark'?'Clair':'Sombre'; } }catch(e){}
pushToInputs();
renderAll();
