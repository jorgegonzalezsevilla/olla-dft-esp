(function () {
'use strict';
const data=JSON.parse(document.getElementById('studio-data').textContent);
const $=id=>document.getElementById(id), ns='http://www.w3.org/2000/svg';
let lang=data.language, page=0, selected=new Set(data.rows.map(r=>r.id)), visible=[], chosen=[], plotted=[], chartOK=false, busy=false;
const defaults={x:'case',y:'',title:'',kind:'points',color:'#0072B2',size:'report',font:'22',min:'',max:''};
let view={...defaults}, saved=data.view;
const t=k=>(data.labels[lang]||data.labels.en)[k]||k;
const finite=x=>typeof x==='number'&&Number.isFinite(x);
const fmt=x=>!finite(x)?'—':(x!==0&&(Math.abs(x)<0.001||Math.abs(x)>=100000)?x.toExponential(4):Number(x.toPrecision(7)).toString());
const key=(metric,unit)=>JSON.stringify([metric,unit]);
const axes=new Map();
for(const row of data.rows)for(const [name,m] of Object.entries(row.metrics))axes.set(key(name,m.unit),{name,unit:m.unit});
const axisLabel=k=>k==='case'?t('recordAxis'):(axes.has(k)?t(axes.get(k).name)+' ['+(axes.get(k).unit||t('unitMissing'))+']':'—');
const value=(row,k)=>{if(k==='case')return row.index;const a=axes.get(k),m=a&&row.metrics[a.name];return m&&m.unit===a.unit&&finite(m.value)?m.value:null;};
const option=(value,text)=>{const n=document.createElement('option');n.value=value;n.textContent=text;return n;};
function labels(){
 document.documentElement.lang=lang;
 document.querySelectorAll('[data-label]').forEach(n=>{n.textContent=t(n.dataset.label);});
 document.querySelectorAll('[data-aria]').forEach(n=>n.setAttribute('aria-label',t(n.dataset.aria)));
 const calc=$('calculation').value;$('calculation').replaceChildren(option('all',t('all')),...[...new Set(data.rows.map(r=>r.calculation))].sort().map(v=>option(v,v||'—')));$('calculation').value=calc||'all';
 $('metric-y').replaceChildren(...[...axes].map(([k])=>option(k,axisLabel(k))));
 $('metric-x').replaceChildren(option('case',t('record')),...[...axes].map(([k])=>option(k,axisLabel(k))));
 $('metric-y').value=view.y;$('metric-x').value=view.x;
 $('chart-title').placeholder=t('defaultTitle');$('language').value=lang;
 $('loaded-count').textContent=data.rows.length+' '+t('loaded')+' / '+data.total_count;
 $('generated').textContent=data.generated;$('version').textContent=data.qekit_version;
 $('project-title').textContent=data.title;
 $('limit-note').hidden=data.rows.length>=data.total_count;
 $('limit-note').textContent=t(data.exported_scope?'subset':'truncated');
 const present=new Set(data.rows.map(r=>r.id)),missing=(data.view?.selected||[]).filter(id=>!present.has(id)).length;
 $('restore-note').hidden=missing===0;$('restore-note').textContent=missing+' '+t('missingSelection');
}
const available=[...axes.keys()].filter(k=>data.rows.some(r=>finite(value(r,k))));
view.y=available.find(k=>axes.get(k).name==='energy_per_atom')||available[0]||[...axes.keys()][0]||'';
if(saved){view={...view,...saved.figure};selected=new Set(saved.selected||[]);}
labels();
if(saved){$('search').value=saved.search||'';$('status').value=saved.status||'converged';$('calculation').value=saved.calculation||'all';$('plot-only').checked=saved.plotOnly!==false;$('scale').value=saved.scale||'2';}
function syncControls(){for(const [id,k] of Object.entries({'metric-x':'x','metric-y':'y','chart-title':'title','kind':'kind','color':'color','size':'size','font':'font','y-min':'min','y-max':'max'}))$(id).value=view[k];}
syncControls();
function node(tag,attrs={},text){const n=document.createElementNS(ns,tag);for(const [k,v] of Object.entries(attrs))n.setAttribute(k,String(v));if(text!==undefined)n.textContent=text;return n;}
function words(text,max){const parts=[];let rest=String(text);while(rest.length>max&&parts.length<1){let at=rest.lastIndexOf(' ',max);if(at<max/3)at=max;parts.push(rest.slice(0,at));rest=rest.slice(at).trim();}parts.push(rest.length>max?rest.slice(0,max-1)+'…':rest);return parts;}
function draw(){
 plotted=[];chartOK=false;$('chart-wrap').replaceChildren();
 const title=view.title.trim()||t('defaultTitle');$('figure-heading').textContent=title;
 const sizes={report:[960,640],slide:[1280,720],square:[960,960]},[w,h]=sizes[view.size]||sizes.report,fs=Number(view.font)||22;
 const eligible=chosen.filter(r=>finite(value(r,view.y))&&finite(value(r,view.x)));
 let error='';
 if(!eligible.length)error=t('noPoints');
 if(eligible.length>2000)error=t('tooMany');
 let lo,hi,xlo,xhi;
 const bar=view.kind==='bars'&&view.x==='case';
 if(!error){
  const ys=eligible.map(r=>value(r,view.y)),xs=eligible.map(r=>value(r,view.x));
  lo=Math.min(...ys);hi=Math.max(...ys);xlo=Math.min(...xs);xhi=Math.max(...xs);
  const pad=Math.max((hi-lo)*.09,Math.abs(hi)*.000001,1e-8);lo-=pad;hi+=pad;
  if(bar){lo=Math.min(0,lo);hi=Math.max(0,hi);}
  if(view.min!=='')lo=Number(view.min);if(view.max!=='')hi=Number(view.max);
  if($('y-min').validity.badInput||$('y-max').validity.badInput||!finite(lo)||!finite(hi)||!finite(hi-lo)||lo>=hi)error=t('badRange');
  if(!error&&bar&&(lo>0||hi<0))error=t('barZero');
  if(xlo===xhi){xlo-=.5;xhi+=.5;}else{const padX=(xhi-xlo)*.04;xlo-=padX;xhi+=padX;}
  if(!finite(xlo)||!finite(xhi)||!finite(xhi-xlo))error=t('badRange');
 }
 const warnings=[];
 if(chosen.some(r=>r.status!=='converged'))warnings.push(t('stateCaution'));
 if(bar&&axes.get(view.y)?.name.startsWith('energy'))warnings.push(t('energyBars'));
 if(error)warnings.push(error);
 else {
  plotted=eligible.filter(r=>value(r,view.y)>=lo&&value(r,view.y)<=hi);
  if(!plotted.length){error=t('noPoints');warnings.push(error);}
 }
 if(chosen.some(r=>Object.values(r.metrics).some(m=>finite(m.uncertainty))))warnings.push(t('uncertaintyNote'));
 const distinct=new Set(chosen.map(r=>r.formula+'|'+r.method_sha256));
 if(chosen.length&&(distinct.size>1||chosen.some(r=>!r.method_known||!r.parameters||!Object.keys(r.parameters).length)))warnings.push(t('mixed'));
 const omitted=chosen.length-plotted.length;
 if(omitted&&!error)warnings.push(omitted+' '+t('omitted')+'. '+t('units'));
 $('plot-warning').textContent=warnings.join(' ');$('point-count').textContent=plotted.length+' '+t('plotted');
 $('point-info').textContent=t('pointHint');
 if(error){const p=document.createElement('p');p.style.padding='45px 24px';p.textContent=error;$('chart-wrap').append(p);return;}
 const svg=node('svg',{xmlns:ns,viewBox:`0 0 ${w} ${h}`,width:w,height:h,role:'img','aria-labelledby':'plot-title plot-description'});
 svg.append(node('title',{id:'plot-title'},title),node('desc',{id:'plot-description'},axisLabel(view.y)+' · '+axisLabel(view.x)+'. '+warnings.join(' ')));
 svg.append(node('rect',{x:0,y:0,width:w,height:h,fill:'#ffffff'}));
 const ink='#203b42',muted='#526a70',left=130,right=w-45,top=120,bottom=h-130;
 function text(x,y,s,attrs={}){svg.append(node('text',{x,y,fill:ink,'font-family':'Arial, sans-serif','font-size':fs,...attrs},s));}
 words(title,Math.floor((w-100)/(fs*.57))).forEach((line,i)=>text(44,42+i*(fs+6),line,{'font-size':fs+3,'font-weight':600}));
 const offsetY=Math.abs((hi+lo)/2)>0&&(hi-lo)/Math.abs((hi+lo)/2)<1e-4?value(eligible[0],view.y):0;
 const offsetX=view.x!=='case'&&Math.abs((xhi+xlo)/2)>0&&(xhi-xlo)/Math.abs((xhi+xlo)/2)<1e-4?value(eligible[0],view.x):0;
 const offsetNote=[offsetY?'Y = '+String(offsetY)+' + ΔY':'',offsetX?'X = '+String(offsetX)+' + ΔX':''].filter(Boolean).join(' · ');
 text(44,99,offsetNote||t('rangeHint')+(distinct.size>1?' · '+t('exploratory'):''),{'font-size':14,fill:muted});
 const px=x=>left+(x-xlo)/(xhi-xlo)*(right-left),py=y=>bottom-(y-lo)/(hi-lo)*(bottom-top);
 for(let i=0;i<=5;i++){
  const y=lo+(hi-lo)*i/5,yy=py(y);
  svg.append(node('line',{x1:left,y1:yy,x2:right,y2:yy,stroke:'#dce5e7','stroke-width':1}));
  text(left-12,yy+5,fmt(y-offsetY),{'text-anchor':'end','font-size':fs-4,fill:muted});
 }
 svg.append(node('line',{x1:left,y1:top,x2:left,y2:bottom,stroke:'#90a8ae'}),node('line',{x1:left,y1:bottom,x2:right,y2:bottom,stroke:'#90a8ae'}));
 if(view.x==='case'){
  const stride=Math.max(1,Math.ceil(eligible.length/10));eligible.forEach((r,i)=>{if(i%stride===0)text(px(r.index),bottom+30,String(r.index),{'text-anchor':'middle','font-size':fs-3,fill:muted});});
 }else for(let i=0;i<=5;i++)text(left+(right-left)*i/5,bottom+30,fmt(xlo+(xhi-xlo)*i/5-offsetX),{'text-anchor':'middle','font-size':fs-4,fill:muted});
 text((left+right)/2,h-67,(offsetX?'Δ ':'')+axisLabel(view.x),{'text-anchor':'middle'});
 const yl=node('text',{transform:`translate(30 ${(top+bottom)/2}) rotate(-90)`,'text-anchor':'middle',fill:ink,'font-family':'Arial, sans-serif','font-size':fs},(offsetY?'Δ ':'')+axisLabel(view.y));svg.append(yl);
 const barW=Math.max(1,Math.min(42,(right-left)/Math.max(eligible.length,1)*.65));
 for(const r of plotted){
  const x=px(value(r,view.x)),y=py(value(r,view.y)),ok=r.status==='converged';
  const description=`${t('record')} ${r.index} · ${r.formula} · ${axisLabel(view.y)}: ${String(value(r,view.y))} · ${t(r.status)} · ${r.id}`;
  const group=node('g',{tabindex:0,role:'img','aria-label':description});group.append(node('title',{},description));
  if(bar)group.append(node('rect',{x:x-barW/2,y:Math.min(y,py(0)),width:barW,height:Math.max(.6,Math.abs(y-py(0))),fill:ok?view.color:'#ffffff',stroke:view.color,'stroke-width':2,'stroke-dasharray':ok?'none':'4 3'}));
  else if(ok)group.append(node('circle',{cx:x,cy:y,r:6,fill:view.color,stroke:'#fff','stroke-width':1}));
  else group.append(node('path',{d:`M${x-6},${y-6}L${x+6},${y+6}M${x+6},${y-6}L${x-6},${y+6}`,stroke:view.color,'stroke-width':3,fill:'none'}));
  group.addEventListener('focus',()=>{$('point-info').textContent=description;});group.addEventListener('pointerenter',()=>{$('point-info').textContent=description;});svg.append(group);
 }
 text(left,h-32,'● '+t('converged')+'    × / ▱ '+t('otherStates'),{'font-size':13,fill:muted});
 text(left,h-12,plotted.some(r=>r.status!=='converged')?t('stateFigure'):plotted.length+' '+t('plotted')+' · '+omitted+' '+t('omitted'),{'font-size':11,fill:muted});
 $('chart-wrap').append(svg);chartOK=true;
}
function table(){
 $('rows').replaceChildren();const pages=Math.max(1,Math.ceil(visible.length/50));page=Math.min(page,pages-1);
 for(const r of visible.slice(page*50,(page+1)*50)){
  const tr=document.createElement('tr');
  function cell(text){const c=document.createElement('td');if(text!==undefined)c.textContent=text;tr.append(c);return c;}
  const check=document.createElement('input');check.type='checkbox';check.checked=selected.has(r.id);check.setAttribute('aria-label',t('use')+' '+t('record')+' '+r.index+' '+r.formula);
  check.addEventListener('change',()=>{check.checked?selected.add(r.id):selected.delete(r.id);render();const again=$('rows').querySelector(`[data-row="${r.index}"]`);if(again)again.focus();});check.dataset.row=r.index;cell().append(check);
  const id=cell(String(r.index));const small=document.createElement('small');small.textContent=r.id.slice(0,10);id.append(small);cell(r.formula||'—');const tag=cell(r.tag||'—');tag.className='tag-cell';tag.title=r.tag;
  const badge=document.createElement('span');badge.className='badge'+(r.status==='converged'?' ok':'');badge.textContent=(r.status==='converged'?'● ':'× ')+t(r.status);cell().append(badge);
  const val=value(r,view.y);const v=cell(fmt(val));v.title=finite(val)?String(val):t('missing');cell(t(r.review));
  const detail=document.createElement('button');detail.textContent=t('details');detail.setAttribute('aria-label',t('details')+' '+t('record')+' '+r.index);detail.addEventListener('click',()=>{$('record-details').hidden=false;$('detail-text').textContent=JSON.stringify(r,null,2);$('record-details').scrollIntoView({block:'nearest'});});cell().append(detail);$('rows').append(tr);
 }
 $('value-heading').textContent=axisLabel(view.y);$('empty-table').hidden=visible.length>0;
 $('page-count').textContent=t('page')+' '+(page+1)+' / '+pages;$('previous').disabled=page===0;$('next').disabled=page>=pages-1;
}
function exportRows(){return $('plot-only').checked?plotted:chosen;}
function scope(){const n=exportRows().length;$('export-scope').textContent=t('exporting')+' '+n+' '+t('records')+'. '+t('scope');for(const id of ['save-csv','save-json','save-html'])$(id).disabled=n===0;$('save-svg').disabled=!chartOK;$('save-png').disabled=!chartOK||busy;}
function render(){
 const q=$('search').value.trim().toLocaleLowerCase(),status=$('status').value,calc=$('calculation').value;
 visible=data.rows.filter(r=>(status==='all'||r.status===status)&&(calc==='all'||r.calculation===calc)&&(!q||[r.formula,r.tag,r.id].join(' ').toLocaleLowerCase().includes(q)));
 chosen=visible.filter(r=>selected.has(r.id));$('selection-status').textContent=visible.length+' '+t('filtered')+' / '+data.rows.length+' '+t('loaded')+' · '+chosen.length+' '+t('selected');
 $('kind').querySelector('[value=bars]').disabled=view.x!=='case';if(view.x!=='case'&&view.kind==='bars'){view.kind='points';$('kind').value='points';}
 draw();table();scope();
}
let pending;
$('search').addEventListener('input',()=>{clearTimeout(pending);pending=setTimeout(()=>{page=0;render();},120);});
for(const id of ['status','calculation'])$(id).addEventListener('change',()=>{page=0;render();});
for(const [id,k] of Object.entries({'metric-x':'x','metric-y':'y','chart-title':'title','kind':'kind','color':'color','size':'size','font':'font','y-min':'min','y-max':'max'}))$(id).addEventListener(id==='chart-title'||id==='y-min'||id==='y-max'?'input':'change',()=>{view[k]=$(id).value;if(id==='metric-y'){view.min='';view.max='';$('y-min').value='';$('y-max').value='';}render();});
$('language').addEventListener('change',()=>{lang=$('language').value;$('download-status').textContent='';labels();render();});
$('reset').addEventListener('click',()=>{view={...defaults,y:view.y};syncControls();render();});
$('zoom-chart').addEventListener('click',()=>{const zoomed=$('chart-wrap').classList.toggle('zoomed');$('zoom-chart').setAttribute('aria-pressed',String(zoomed));});
$('plot-only').addEventListener('change',scope);
$('select-all').addEventListener('click',()=>{visible.forEach(r=>selected.add(r.id));render();});$('clear').addEventListener('click',()=>{selected.clear();render();});
$('previous').addEventListener('click',()=>{page--;table();});$('next').addEventListener('click',()=>{page++;table();});
document.querySelector('.export-buttons').addEventListener('click',()=>{clearTimeout(pending);render();},true);
function state(){return {figure:{...view},search:$('search').value,status:$('status').value,calculation:$('calculation').value,plotOnly:$('plot-only').checked,scale:$('scale').value,selected:exportRows().map(r=>r.id)};}
function name(ext){return 'olla-results.'+ext;}
function download(blob,filename){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=filename;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),2000);$('download-status').textContent=t('saved')+': '+filename;}
function payload(){return {schema_version:1,order:data.order,qekit_version:data.qekit_version,generated:data.generated,exported:new Date().toISOString(),source_count:data.total_count,loaded_count:data.rows.length,count:exportRows().length,view:state(),results:exportRows()};}
function csvCell(v){if(v===null||v===undefined)return '';let s=String(v);if(typeof v!=='number'&&/^[\s\u0000-\u001f]*[=+@-]/.test(s))s="'"+s;return '"'+s.replace(/"/g,'""')+'"';}
$('save-json').addEventListener('click',()=>download(new Blob([JSON.stringify(payload(),null,2)],{type:'application/json;charset=utf-8'}),name('json')));
$('save-csv').addEventListener('click',()=>{const rows=exportRows(),metrics=[...new Set(rows.flatMap(r=>Object.keys(r.metrics)))].sort(),base=['id','index','formula','calculation','status','review','tag','source_sha256','method_sha256'],header=[...base,...metrics.flatMap(k=>[k+'.value',k+'.unit',k+'.reason',k+'.uncertainty'])];const lines=[header.map(csvCell).join(',')];for(const r of rows)lines.push([...base.map(k=>r[k]),...metrics.flatMap(k=>{const m=r.metrics[k]||{};return [m.value,m.unit,m.reason,m.uncertainty];})].map(csvCell).join(','));download(new Blob(['\ufeff'+lines.join('\r\n')+'\r\n'],{type:'text/csv;charset=utf-8'}),name('csv'));});
function svgText(){const clone=$('chart-wrap').querySelector('svg').cloneNode(true);clone.querySelectorAll('[tabindex]').forEach(n=>n.removeAttribute('tabindex'));const metadata=node('metadata',{},JSON.stringify({schema_version:1,qekit_version:data.qekit_version,generated:data.generated,exported:new Date().toISOString(),order:data.order,view:state(),ids:plotted.map(r=>r.id),units:{x:axisLabel(view.x),y:axisLabel(view.y)},method_hashes:[...new Set(plotted.map(r=>r.method_sha256))]}));clone.append(metadata);return new XMLSerializer().serializeToString(clone);}
$('save-svg').addEventListener('click',()=>{if(chartOK)download(new Blob([svgText()],{type:'image/svg+xml;charset=utf-8'}),name('svg'));});
$('save-png').addEventListener('click',async()=>{
 if(!chartOK||busy)return;busy=true;scope();$('download-status').textContent=t('pngBusy');let url;
 try{const svg=$('chart-wrap').querySelector('svg'),scale=Number($('scale').value),w=Number(svg.getAttribute('width'))*scale,h=Number(svg.getAttribute('height'))*scale;
  if(![1,2].includes(scale)||w*h>6000000)throw Error('Image limit');
  url=URL.createObjectURL(new Blob([svgText()],{type:'image/svg+xml;charset=utf-8'}));const img=new Image();await new Promise((resolve,reject)=>{img.onload=resolve;img.onerror=reject;img.src=url;});
  const canvas=document.createElement('canvas');canvas.width=w;canvas.height=h;const ctx=canvas.getContext('2d');ctx.drawImage(img,0,0,w,h);
  const blob=await new Promise(resolve=>canvas.toBlob(resolve,'image/png'));if(!blob)throw Error('PNG unavailable');download(blob,name('png'));
 }catch(_){$('download-status').textContent=t('failed');}finally{if(url)URL.revokeObjectURL(url);busy=false;scope();}
});
$('save-html').addEventListener('click',()=>{
 const snapshot={...data,exported_scope:true,language:lang,rows:exportRows(),view:state()};const clone=document.documentElement.cloneNode(true);
 clone.querySelector('#studio-data').textContent=JSON.stringify(snapshot).replace(/</g,'\\u003c').replace(/>/g,'\\u003e').replace(/&/g,'\\u0026');
 clone.querySelector('#point-info').textContent='';clone.querySelector('#rows').replaceChildren();clone.querySelector('#chart-wrap').replaceChildren();clone.querySelector('#detail-text').textContent='';clone.querySelector('#record-details').hidden=true;clone.querySelector('#download-status').textContent='';
 download(new Blob(['<!doctype html>\n'+clone.outerHTML],{type:'text/html;charset=utf-8'}),name('html'));
});
render();
})();
