var __defProp = Object.defineProperty;
var __name = (target, value) => __defProp(target, "name", { value, configurable: true });
import wk from "./node-worker";
const u8 = Uint8Array, u16 = Uint16Array, i32 = Int32Array;
const fleb = new u8([
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  1,
  1,
  1,
  1,
  2,
  2,
  2,
  2,
  3,
  3,
  3,
  3,
  4,
  4,
  4,
  4,
  5,
  5,
  5,
  5,
  0,
  /* unused */
  0,
  0,
  /* impossible */
  0
]);
const fdeb = new u8([
  0,
  0,
  0,
  0,
  1,
  1,
  2,
  2,
  3,
  3,
  4,
  4,
  5,
  5,
  6,
  6,
  7,
  7,
  8,
  8,
  9,
  9,
  10,
  10,
  11,
  11,
  12,
  12,
  13,
  13,
  /* unused */
  0,
  0
]);
const clim = new u8([16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15]);
const freb = /* @__PURE__ */ __name((eb, start) => {
  const b = new u16(31);
  for (let i = 0; i < 31; ++i) {
    b[i] = start += 1 << eb[i - 1];
  }
  const r = new i32(b[30]);
  for (let i = 1; i < 30; ++i) {
    for (let j = b[i]; j < b[i + 1]; ++j) {
      r[j] = j - b[i] << 5 | i;
    }
  }
  return { b, r };
}, "freb");
const { b: fl, r: revfl } = freb(fleb, 2);
fl[28] = 258, revfl[258] = 28;
const { b: fd, r: revfd } = freb(fdeb, 0);
const rev = new u16(32768);
for (let i = 0; i < 32768; ++i) {
  let x = (i & 43690) >> 1 | (i & 21845) << 1;
  x = (x & 52428) >> 2 | (x & 13107) << 2;
  x = (x & 61680) >> 4 | (x & 3855) << 4;
  rev[i] = ((x & 65280) >> 8 | (x & 255) << 8) >> 1;
}
const hMap = /* @__PURE__ */ __name((cd, mb, r) => {
  const s = cd.length;
  let i = 0;
  const l = new u16(mb);
  for (; i < s; ++i) {
    if (cd[i]) ++l[cd[i] - 1];
  }
  const le = new u16(mb);
  for (i = 1; i < mb; ++i) {
    le[i] = le[i - 1] + l[i - 1] << 1;
  }
  let co;
  if (r) {
    co = new u16(1 << mb);
    const rvb = 15 - mb;
    for (i = 0; i < s; ++i) {
      if (cd[i]) {
        const sv = i << 4 | cd[i];
        const r2 = mb - cd[i];
        let v = le[cd[i] - 1]++ << r2;
        for (const m = v | (1 << r2) - 1; v <= m; ++v) {
          co[rev[v] >> rvb] = sv;
        }
      }
    }
  } else {
    co = new u16(s);
    for (i = 0; i < s; ++i) {
      if (cd[i]) {
        co[i] = rev[le[cd[i] - 1]++] >> 15 - cd[i];
      }
    }
  }
  return co;
}, "hMap");
const flt = new u8(288);
for (let i = 0; i < 144; ++i) flt[i] = 8;
for (let i = 144; i < 256; ++i) flt[i] = 9;
for (let i = 256; i < 280; ++i) flt[i] = 7;
for (let i = 280; i < 288; ++i) flt[i] = 8;
const fdt = new u8(32);
for (let i = 0; i < 32; ++i) fdt[i] = 5;
const flm = /* @__PURE__ */ hMap(flt, 9, 0), flrm = /* @__PURE__ */ hMap(flt, 9, 1);
const fdm = /* @__PURE__ */ hMap(fdt, 5, 0), fdrm = /* @__PURE__ */ hMap(fdt, 5, 1);
const max = /* @__PURE__ */ __name((a) => {
  let m = a[0];
  for (let i = 1; i < a.length; ++i) {
    if (a[i] > m) m = a[i];
  }
  return m;
}, "max");
const bits = /* @__PURE__ */ __name((d, p, m) => {
  const o = p / 8 | 0;
  return (d[o] | d[o + 1] << 8) >> (p & 7) & m;
}, "bits");
const bits16 = /* @__PURE__ */ __name((d, p) => {
  const o = p / 8 | 0;
  return (d[o] | d[o + 1] << 8 | d[o + 2] << 16) >> (p & 7);
}, "bits16");
const shft = /* @__PURE__ */ __name((p) => (p + 7) / 8 | 0, "shft");
const slc = /* @__PURE__ */ __name((v, s, e) => {
  if (s == null || s < 0) s = 0;
  if (e == null || e > v.length) e = v.length;
  return new u8(v.subarray(s, e));
}, "slc");
const FlateErrorCode = {
  UnexpectedEOF: 0,
  InvalidBlockType: 1,
  InvalidLengthLiteral: 2,
  InvalidDistance: 3,
  StreamFinished: 4,
  NoStreamHandler: 5,
  InvalidHeader: 6,
  NoCallback: 7,
  InvalidUTF8: 8,
  ExtraFieldTooLong: 9,
  InvalidDate: 10,
  FilenameTooLong: 11,
  StreamFinishing: 12,
  InvalidZipData: 13,
  UnknownCompressionMethod: 14
};
const ec = [
  "unexpected EOF",
  "invalid block type",
  "invalid length/literal",
  "invalid distance",
  "stream finished",
  "no stream handler",
  ,
  // determined by compression function
  "no callback",
  "invalid UTF-8 data",
  "extra field too long",
  "date not in range 1980-2099",
  "filename too long",
  "stream finishing",
  "invalid zip data"
  // determined by unknown compression method
];
;
const err = /* @__PURE__ */ __name((ind, msg, nt) => {
  const e = new Error(msg || ec[ind]);
  e.code = ind;
  if (Error.captureStackTrace) Error.captureStackTrace(e, err);
  if (!nt) throw e;
  return e;
}, "err");
const inflt = /* @__PURE__ */ __name((dat, st, buf, dict) => {
  const sl = dat.length, dl = dict ? dict.length : 0;
  if (!sl || st.f && !st.l) return buf || new u8(0);
  const noBuf = !buf;
  const resize = noBuf || st.i != 2;
  const noSt = st.i;
  if (noBuf) buf = new u8(sl * 3);
  const cbuf = /* @__PURE__ */ __name((l) => {
    let bl = buf.length;
    if (l > bl) {
      const nbuf = new u8(Math.max(bl * 2, l));
      nbuf.set(buf);
      buf = nbuf;
    }
  }, "cbuf");
  let final = st.f || 0, pos = st.p || 0, bt = st.b || 0, lm = st.l, dm = st.d, lbt = st.m, dbt = st.n;
  const tbts = sl * 8;
  do {
    if (!lm) {
      final = bits(dat, pos, 1);
      const type = bits(dat, pos + 1, 3);
      pos += 3;
      if (!type) {
        const s = shft(pos) + 4, l = dat[s - 4] | dat[s - 3] << 8, t = s + l;
        if (t > sl) {
          if (noSt) err(0);
          break;
        }
        if (resize) cbuf(bt + l);
        buf.set(dat.subarray(s, t), bt);
        st.b = bt += l, st.p = pos = t * 8, st.f = final;
        continue;
      } else if (type == 1) lm = flrm, dm = fdrm, lbt = 9, dbt = 5;
      else if (type == 2) {
        const hLit = bits(dat, pos, 31) + 257, hcLen = bits(dat, pos + 10, 15) + 4;
        const tl = hLit + bits(dat, pos + 5, 31) + 1;
        pos += 14;
        const ldt = new u8(tl);
        const clt = new u8(19);
        for (let i = 0; i < hcLen; ++i) {
          clt[clim[i]] = bits(dat, pos + i * 3, 7);
        }
        pos += hcLen * 3;
        const clb = max(clt), clbmsk = (1 << clb) - 1;
        const clm = hMap(clt, clb, 1);
        for (let i = 0; i < tl; ) {
          const r = clm[bits(dat, pos, clbmsk)];
          pos += r & 15;
          const s = r >> 4;
          if (s < 16) {
            ldt[i++] = s;
          } else {
            let c = 0, n = 0;
            if (s == 16) n = 3 + bits(dat, pos, 3), pos += 2, c = ldt[i - 1];
            else if (s == 17) n = 3 + bits(dat, pos, 7), pos += 3;
            else if (s == 18) n = 11 + bits(dat, pos, 127), pos += 7;
            while (n--) ldt[i++] = c;
          }
        }
        const lt = ldt.subarray(0, hLit), dt = ldt.subarray(hLit);
        lbt = max(lt);
        dbt = max(dt);
        lm = hMap(lt, lbt, 1);
        dm = hMap(dt, dbt, 1);
      } else err(1);
      if (pos > tbts) {
        if (noSt) err(0);
        break;
      }
    }
    if (resize) cbuf(bt + 131072);
    const lms = (1 << lbt) - 1, dms = (1 << dbt) - 1;
    let lpos = pos;
    for (; ; lpos = pos) {
      const c = lm[bits16(dat, pos) & lms], sym = c >> 4;
      pos += c & 15;
      if (pos > tbts) {
        if (noSt) err(0);
        break;
      }
      if (!c) err(2);
      if (sym < 256) buf[bt++] = sym;
      else if (sym == 256) {
        lpos = pos, lm = null;
        break;
      } else {
        let add = sym - 254;
        if (sym > 264) {
          const i = sym - 257, b = fleb[i];
          add = bits(dat, pos, (1 << b) - 1) + fl[i];
          pos += b;
        }
        const d = dm[bits16(dat, pos) & dms], dsym = d >> 4;
        if (!d) err(3);
        pos += d & 15;
        let dt = fd[dsym];
        if (dsym > 3) {
          const b = fdeb[dsym];
          dt += bits16(dat, pos) & (1 << b) - 1, pos += b;
        }
        if (pos > tbts) {
          if (noSt) err(0);
          break;
        }
        if (resize) cbuf(bt + 131072);
        const end = bt + add;
        if (bt < dt) {
          const shift = dl - dt, dend = Math.min(dt, end);
          if (shift + bt < 0) err(3);
          for (; bt < dend; ++bt) buf[bt] = dict[shift + bt];
        }
        for (; bt < end; ++bt) buf[bt] = buf[bt - dt];
      }
    }
    st.l = lm, st.p = lpos, st.b = bt, st.f = final;
    if (lm) final = 1, st.m = lbt, st.d = dm, st.n = dbt;
  } while (!final);
  return bt != buf.length && noBuf ? slc(buf, 0, bt) : buf.subarray(0, bt);
}, "inflt");
const wbits = /* @__PURE__ */ __name((d, p, v) => {
  v <<= p & 7;
  const o = p / 8 | 0;
  d[o] |= v;
  d[o + 1] |= v >> 8;
}, "wbits");
const wbits16 = /* @__PURE__ */ __name((d, p, v) => {
  v <<= p & 7;
  const o = p / 8 | 0;
  d[o] |= v;
  d[o + 1] |= v >> 8;
  d[o + 2] |= v >> 16;
}, "wbits16");
const hTree = /* @__PURE__ */ __name((d, mb) => {
  const t = [];
  for (let i = 0; i < d.length; ++i) {
    if (d[i]) t.push({ s: i, f: d[i] });
  }
  const s = t.length;
  const t2 = t.slice();
  if (!s) return { t: et, l: 0 };
  if (s == 1) {
    const v = new u8(t[0].s + 1);
    v[t[0].s] = 1;
    return { t: v, l: 1 };
  }
  t.sort((a, b) => a.f - b.f);
  t.push({ s: -1, f: 25001 });
  let l = t[0], r = t[1], i0 = 0, i1 = 1, i2 = 2;
  t[0] = { s: -1, f: l.f + r.f, l, r };
  while (i1 != s - 1) {
    l = t[t[i0].f < t[i2].f ? i0++ : i2++];
    r = t[i0 != i1 && t[i0].f < t[i2].f ? i0++ : i2++];
    t[i1++] = { s: -1, f: l.f + r.f, l, r };
  }
  let maxSym = t2[0].s;
  for (let i = 1; i < s; ++i) {
    if (t2[i].s > maxSym) maxSym = t2[i].s;
  }
  const tr = new u16(maxSym + 1);
  let mbt = ln(t[i1 - 1], tr, 0);
  if (mbt > mb) {
    let i = 0, dt = 0;
    const lft = mbt - mb, cst = 1 << lft;
    t2.sort((a, b) => tr[b.s] - tr[a.s] || a.f - b.f);
    for (; i < s; ++i) {
      const i22 = t2[i].s;
      if (tr[i22] > mb) {
        dt += cst - (1 << mbt - tr[i22]);
        tr[i22] = mb;
      } else break;
    }
    dt >>= lft;
    while (dt > 0) {
      const i22 = t2[i].s;
      if (tr[i22] < mb) dt -= 1 << mb - tr[i22]++ - 1;
      else ++i;
    }
    for (; i >= 0 && dt; --i) {
      const i22 = t2[i].s;
      if (tr[i22] == mb) {
        --tr[i22];
        ++dt;
      }
    }
    mbt = mb;
  }
  return { t: new u8(tr), l: mbt };
}, "hTree");
const ln = /* @__PURE__ */ __name((n, l, d) => {
  return n.s == -1 ? Math.max(ln(n.l, l, d + 1), ln(n.r, l, d + 1)) : l[n.s] = d;
}, "ln");
const lc = /* @__PURE__ */ __name((c) => {
  let s = c.length;
  while (s && !c[--s]) ;
  const cl = new u16(++s);
  let cli = 0, cln = c[0], cls = 1;
  const w = /* @__PURE__ */ __name((v) => {
    cl[cli++] = v;
  }, "w");
  for (let i = 1; i <= s; ++i) {
    if (c[i] == cln && i != s)
      ++cls;
    else {
      if (!cln && cls > 2) {
        for (; cls > 138; cls -= 138) w(32754);
        if (cls > 2) {
          w(cls > 10 ? cls - 11 << 5 | 28690 : cls - 3 << 5 | 12305);
          cls = 0;
        }
      } else if (cls > 3) {
        w(cln), --cls;
        for (; cls > 6; cls -= 6) w(8304);
        if (cls > 2) w(cls - 3 << 5 | 8208), cls = 0;
      }
      while (cls--) w(cln);
      cls = 1;
      cln = c[i];
    }
  }
  return { c: cl.subarray(0, cli), n: s };
}, "lc");
const clen = /* @__PURE__ */ __name((cf, cl) => {
  let l = 0;
  for (let i = 0; i < cl.length; ++i) l += cf[i] * cl[i];
  return l;
}, "clen");
const wfblk = /* @__PURE__ */ __name((out, pos, dat) => {
  const s = dat.length;
  const o = shft(pos + 2);
  out[o] = s & 255;
  out[o + 1] = s >> 8;
  out[o + 2] = out[o] ^ 255;
  out[o + 3] = out[o + 1] ^ 255;
  for (let i = 0; i < s; ++i) out[o + i + 4] = dat[i];
  return (o + 4 + s) * 8;
}, "wfblk");
const wblk = /* @__PURE__ */ __name((dat, out, final, syms, lf, df, eb, li, bs, bl, p) => {
  wbits(out, p++, final);
  ++lf[256];
  const { t: dlt, l: mlb } = hTree(lf, 15);
  const { t: ddt, l: mdb } = hTree(df, 15);
  const { c: lclt, n: nlc } = lc(dlt);
  const { c: lcdt, n: ndc } = lc(ddt);
  const lcfreq = new u16(19);
  for (let i = 0; i < lclt.length; ++i) ++lcfreq[lclt[i] & 31];
  for (let i = 0; i < lcdt.length; ++i) ++lcfreq[lcdt[i] & 31];
  const { t: lct, l: mlcb } = hTree(lcfreq, 7);
  let nlcc = 19;
  for (; nlcc > 4 && !lct[clim[nlcc - 1]]; --nlcc) ;
  const flen = bl + 5 << 3;
  const ftlen = clen(lf, flt) + clen(df, fdt) + eb;
  const dtlen = clen(lf, dlt) + clen(df, ddt) + eb + 14 + 3 * nlcc + clen(lcfreq, lct) + 2 * lcfreq[16] + 3 * lcfreq[17] + 7 * lcfreq[18];
  if (bs >= 0 && flen <= ftlen && flen <= dtlen) return wfblk(out, p, dat.subarray(bs, bs + bl));
  let lm, ll, dm, dl;
  wbits(out, p, 1 + (dtlen < ftlen)), p += 2;
  if (dtlen < ftlen) {
    lm = hMap(dlt, mlb, 0), ll = dlt, dm = hMap(ddt, mdb, 0), dl = ddt;
    const llm = hMap(lct, mlcb, 0);
    wbits(out, p, nlc - 257);
    wbits(out, p + 5, ndc - 1);
    wbits(out, p + 10, nlcc - 4);
    p += 14;
    for (let i = 0; i < nlcc; ++i) wbits(out, p + 3 * i, lct[clim[i]]);
    p += 3 * nlcc;
    const lcts = [lclt, lcdt];
    for (let it = 0; it < 2; ++it) {
      const clct = lcts[it];
      for (let i = 0; i < clct.length; ++i) {
        const len = clct[i] & 31;
        wbits(out, p, llm[len]), p += lct[len];
        if (len > 15) wbits(out, p, clct[i] >> 5 & 127), p += clct[i] >> 12;
      }
    }
  } else {
    lm = flm, ll = flt, dm = fdm, dl = fdt;
  }
  for (let i = 0; i < li; ++i) {
    const sym = syms[i];
    if (sym > 255) {
      const len = sym >> 18 & 31;
      wbits16(out, p, lm[len + 257]), p += ll[len + 257];
      if (len > 7) wbits(out, p, sym >> 23 & 31), p += fleb[len];
      const dst = sym & 31;
      wbits16(out, p, dm[dst]), p += dl[dst];
      if (dst > 3) wbits16(out, p, sym >> 5 & 8191), p += fdeb[dst];
    } else {
      wbits16(out, p, lm[sym]), p += ll[sym];
    }
  }
  wbits16(out, p, lm[256]);
  return p + ll[256];
}, "wblk");
const deo = /* @__PURE__ */ new i32([65540, 131080, 131088, 131104, 262176, 1048704, 1048832, 2114560, 2117632]);
const et = /* @__PURE__ */ new u8(0);
const dflt = /* @__PURE__ */ __name((dat, lvl, plvl, pre, post, st) => {
  const s = st.z || dat.length;
  const o = new u8(pre + s + 5 * (1 + Math.ceil(s / 7e3)) + post);
  const w = o.subarray(pre, o.length - post);
  const lst = st.l;
  let pos = (st.r || 0) & 7;
  if (lvl) {
    if (pos) w[0] = st.r >> 3;
    const opt = deo[lvl - 1];
    const n = opt >> 13, c = opt & 8191;
    const msk = (1 << plvl) - 1;
    const prev = st.p || new u16(32768), head = st.h || new u16(msk + 1);
    const bs1 = Math.ceil(plvl / 3), bs2 = 2 * bs1;
    const hsh = /* @__PURE__ */ __name((i2) => (dat[i2] ^ dat[i2 + 1] << bs1 ^ dat[i2 + 2] << bs2) & msk, "hsh");
    const syms = new i32(25e3);
    const lf = new u16(288), df = new u16(32);
    let lc2 = 0, eb = 0, i = st.i || 0, li = 0, wi = st.w || 0, bs = 0;
    for (; i + 2 < s; ++i) {
      const hv = hsh(i);
      let imod = i & 32767, pimod = head[hv];
      prev[imod] = pimod;
      head[hv] = imod;
      if (wi <= i) {
        const rem = s - i;
        if ((lc2 > 7e3 || li > 24576) && (rem > 423 || !lst)) {
          pos = wblk(dat, w, 0, syms, lf, df, eb, li, bs, i - bs, pos);
          li = lc2 = eb = 0, bs = i;
          for (let j = 0; j < 286; ++j) lf[j] = 0;
          for (let j = 0; j < 30; ++j) df[j] = 0;
        }
        let l = 2, d = 0, ch2 = c, dif = imod - pimod & 32767;
        if (rem > 2 && hv == hsh(i - dif)) {
          const maxn = Math.min(n, rem) - 1;
          const maxd = Math.min(32767, i);
          const ml = Math.min(258, rem);
          while (dif <= maxd && --ch2 && imod != pimod) {
            if (dat[i + l] == dat[i + l - dif]) {
              let nl = 0;
              for (; nl < ml && dat[i + nl] == dat[i + nl - dif]; ++nl) ;
              if (nl > l) {
                l = nl, d = dif;
                if (nl > maxn) break;
                const mmd = Math.min(dif, nl - 2);
                let md = 0;
                for (let j = 0; j < mmd; ++j) {
                  const ti = i - dif + j & 32767;
                  const pti = prev[ti];
                  const cd = ti - pti & 32767;
                  if (cd > md) md = cd, pimod = ti;
                }
              }
            }
            imod = pimod, pimod = prev[imod];
            dif += imod - pimod & 32767;
          }
        }
        if (d) {
          syms[li++] = 268435456 | revfl[l] << 18 | revfd[d];
          const lin = revfl[l] & 31, din = revfd[d] & 31;
          eb += fleb[lin] + fdeb[din];
          ++lf[257 + lin];
          ++df[din];
          wi = i + l;
          ++lc2;
        } else {
          syms[li++] = dat[i];
          ++lf[dat[i]];
        }
      }
    }
    for (i = Math.max(i, wi); i < s; ++i) {
      syms[li++] = dat[i];
      ++lf[dat[i]];
    }
    pos = wblk(dat, w, lst, syms, lf, df, eb, li, bs, i - bs, pos);
    if (!lst) {
      st.r = pos & 7 | w[pos / 8 | 0] << 3;
      pos -= 7;
      st.h = head, st.p = prev, st.i = i, st.w = wi;
    }
  } else {
    for (let i = st.w || 0; i < s + lst; i += 65535) {
      let e = i + 65535;
      if (e >= s) {
        w[pos / 8 | 0] = lst;
        e = s;
      }
      pos = wfblk(w, pos + 1, dat.subarray(i, e));
    }
    st.i = s;
  }
  return slc(o, 0, pre + shft(pos) + post);
}, "dflt");
const crct = /* @__PURE__ */ (() => {
  const t = new Int32Array(256);
  for (let i = 0; i < 256; ++i) {
    let c = i, k = 9;
    while (--k) c = (c & 1 && -306674912) ^ c >>> 1;
    t[i] = c;
  }
  return t;
})();
const crc = /* @__PURE__ */ __name(() => {
  let c = -1;
  return {
    p(d) {
      let cr = c;
      for (let i = 0; i < d.length; ++i) cr = crct[cr & 255 ^ d[i]] ^ cr >>> 8;
      c = cr;
    },
    d() {
      return ~c;
    }
  };
}, "crc");
const adler = /* @__PURE__ */ __name(() => {
  let a = 1, b = 0;
  return {
    p(d) {
      let n = a, m = b;
      const l = d.length | 0;
      for (let i = 0; i != l; ) {
        const e = Math.min(i + 2655, l);
        for (; i < e; ++i) m += n += d[i];
        n = (n & 65535) + 15 * (n >> 16), m = (m & 65535) + 15 * (m >> 16);
      }
      a = n, b = m;
    },
    d() {
      a %= 65521, b %= 65521;
      return (a & 255) << 24 | (a & 65280) << 8 | (b & 255) << 8 | b >> 8;
    }
  };
}, "adler");
;
const dopt = /* @__PURE__ */ __name((dat, opt, pre, post, st) => {
  if (!st) {
    st = { l: 1 };
    if (opt.dictionary) {
      const dict = opt.dictionary.subarray(-32768);
      const newDat = new u8(dict.length + dat.length);
      newDat.set(dict);
      newDat.set(dat, dict.length);
      dat = newDat;
      st.w = dict.length;
    }
  }
  return dflt(dat, opt.level == null ? 6 : opt.level, opt.mem == null ? st.l ? Math.ceil(Math.max(8, Math.min(13, Math.log(dat.length))) * 1.5) : 20 : 12 + opt.mem, pre, post, st);
}, "dopt");
const mrg = /* @__PURE__ */ __name((a, b) => {
  const o = {};
  for (const k in a) o[k] = a[k];
  for (const k in b) o[k] = b[k];
  return o;
}, "mrg");
const wcln = /* @__PURE__ */ __name((fn, fnStr, td2) => {
  const dt = fn();
  const st = fn.toString();
  const ks = st.slice(st.indexOf("[") + 1, st.lastIndexOf("]")).replace(/\s+/g, "").split(",");
  for (let i = 0; i < dt.length; ++i) {
    let v = dt[i], k = ks[i];
    if (typeof v == "function") {
      fnStr += ";" + k + "=";
      const st2 = v.toString();
      if (v.prototype) {
        if (st2.indexOf("[native code]") != -1) {
          const spInd = st2.indexOf(" ", 8) + 1;
          fnStr += st2.slice(spInd, st2.indexOf("(", spInd));
        } else {
          fnStr += st2;
          for (const t in v.prototype) fnStr += ";" + k + ".prototype." + t + "=" + v.prototype[t].toString();
        }
      } else fnStr += st2;
    } else td2[k] = v;
  }
  return fnStr;
}, "wcln");
const ch = [];
const cbfs = /* @__PURE__ */ __name((v) => {
  const tl = [];
  for (const k in v) {
    if (v[k].buffer) {
      tl.push((v[k] = new v[k].constructor(v[k])).buffer);
    }
  }
  return tl;
}, "cbfs");
const wrkr = /* @__PURE__ */ __name((fns, init, id, cb) => {
  if (!ch[id]) {
    let fnStr = "", td3 = {}, m = fns.length - 1;
    for (let i = 0; i < m; ++i)
      fnStr = wcln(fns[i], fnStr, td3);
    ch[id] = { c: wcln(fns[m], fnStr, td3), e: td3 };
  }
  const td2 = mrg({}, ch[id].e);
  return wk(ch[id].c + ";onmessage=function(e){for(var k in e.data)self[k]=e.data[k];onmessage=" + init.toString() + "}", id, td2, cbfs(td2), cb);
}, "wrkr");
const bInflt = /* @__PURE__ */ __name(() => [u8, u16, i32, fleb, fdeb, clim, fl, fd, flrm, fdrm, rev, ec, hMap, max, bits, bits16, shft, slc, err, inflt, inflateSync, pbf, gopt], "bInflt");
const bDflt = /* @__PURE__ */ __name(() => [u8, u16, i32, fleb, fdeb, clim, revfl, revfd, flm, flt, fdm, fdt, rev, deo, et, hMap, wbits, wbits16, hTree, ln, lc, clen, wfblk, wblk, shft, slc, dflt, dopt, deflateSync, pbf], "bDflt");
const gze = /* @__PURE__ */ __name(() => [gzh, gzhl, wbytes, crc, crct], "gze");
const guze = /* @__PURE__ */ __name(() => [gzs, gzl], "guze");
const zle = /* @__PURE__ */ __name(() => [zlh, wbytes, adler], "zle");
const zule = /* @__PURE__ */ __name(() => [zls], "zule");
const pbf = /* @__PURE__ */ __name((msg) => postMessage(msg, [msg.buffer]), "pbf");
const gopt = /* @__PURE__ */ __name((o) => o && {
  out: o.size && new u8(o.size),
  dictionary: o.dictionary
}, "gopt");
const cbify = /* @__PURE__ */ __name((dat, opts, fns, init, id, cb) => {
  const w = wrkr(
    fns,
    init,
    id,
    (err2, dat2) => {
      w.terminate();
      cb(err2, dat2);
    }
  );
  w.postMessage([dat, opts], opts.consume ? [dat.buffer] : []);
  return () => {
    w.terminate();
  };
}, "cbify");
const astrm = /* @__PURE__ */ __name((strm) => {
  strm.ondata = (dat, final) => postMessage([dat, final], [dat.buffer]);
  return (ev) => {
    if (ev.data.length) {
      strm.push(ev.data[0], ev.data[1]);
      postMessage([ev.data[0].length]);
    } else strm.flush();
  };
}, "astrm");
const astrmify = /* @__PURE__ */ __name((fns, strm, opts, init, id, flush, ext) => {
  let t;
  const w = wrkr(
    fns,
    init,
    id,
    (err2, dat) => {
      if (err2) w.terminate(), strm.ondata.call(strm, err2);
      else if (!Array.isArray(dat)) ext(dat);
      else if (dat.length == 1) {
        strm.queuedSize -= dat[0];
        if (strm.ondrain) strm.ondrain(dat[0]);
      } else {
        if (dat[1]) w.terminate();
        strm.ondata.call(strm, err2, dat[0], dat[1]);
      }
    }
  );
  w.postMessage(opts);
  strm.queuedSize = 0;
  strm.push = (d, f) => {
    if (!strm.ondata) err(5);
    if (t) strm.ondata(err(4, 0, 1), null, !!f);
    strm.queuedSize += d.length;
    w.postMessage([d, t = f], [d.buffer]);
  };
  strm.terminate = () => {
    w.terminate();
  };
  if (flush) {
    strm.flush = () => {
      w.postMessage([]);
    };
  }
}, "astrmify");
const b2 = /* @__PURE__ */ __name((d, b) => d[b] | d[b + 1] << 8, "b2");
const b4 = /* @__PURE__ */ __name((d, b) => (d[b] | d[b + 1] << 8 | d[b + 2] << 16 | d[b + 3] << 24) >>> 0, "b4");
const b8 = /* @__PURE__ */ __name((d, b) => b4(d, b) + b4(d, b + 4) * 4294967296, "b8");
const wbytes = /* @__PURE__ */ __name((d, b, v) => {
  for (; v; ++b) d[b] = v, v >>>= 8;
}, "wbytes");
const gzh = /* @__PURE__ */ __name((c, o) => {
  const fn = o.filename;
  c[0] = 31, c[1] = 139, c[2] = 8, c[8] = o.level < 2 ? 4 : o.level == 9 ? 2 : 0, c[9] = 3;
  if (o.mtime != 0) wbytes(c, 4, Math.floor(new Date(o.mtime || Date.now()) / 1e3));
  if (fn) {
    c[3] = 8;
    for (let i = 0; i <= fn.length; ++i) c[i + 10] = fn.charCodeAt(i);
  }
}, "gzh");
const gzs = /* @__PURE__ */ __name((d) => {
  if (d[0] != 31 || d[1] != 139 || d[2] != 8) err(6, "invalid gzip data");
  const flg = d[3];
  let st = 10;
  if (flg & 4) st += (d[10] | d[11] << 8) + 2;
  for (let zs = (flg >> 3 & 1) + (flg >> 4 & 1); zs > 0; zs -= !d[st++]) ;
  return st + (flg & 2);
}, "gzs");
const gzl = /* @__PURE__ */ __name((d) => {
  const l = d.length;
  return (d[l - 4] | d[l - 3] << 8 | d[l - 2] << 16 | d[l - 1] << 24) >>> 0;
}, "gzl");
const gzhl = /* @__PURE__ */ __name((o) => 10 + (o.filename ? o.filename.length + 1 : 0), "gzhl");
const zlh = /* @__PURE__ */ __name((c, o) => {
  const lv = o.level, fl2 = lv == 0 ? 0 : lv < 6 ? 1 : lv == 9 ? 3 : 2;
  c[0] = 120, c[1] = fl2 << 6 | (o.dictionary && 32);
  c[1] |= 31 - (c[0] << 8 | c[1]) % 31;
  if (o.dictionary) {
    const h = adler();
    h.p(o.dictionary);
    wbytes(c, 2, h.d());
  }
}, "zlh");
const zls = /* @__PURE__ */ __name((d, dict) => {
  if ((d[0] & 15) != 8 || d[0] >> 4 > 7 || (d[0] << 8 | d[1]) % 31) err(6, "invalid zlib data");
  if ((d[1] >> 5 & 1) == +!dict) err(6, "invalid zlib data: " + (d[1] & 32 ? "need" : "unexpected") + " dictionary");
  return (d[1] >> 3 & 4) + 2;
}, "zls");
function StrmOpt(opts, cb) {
  if (typeof opts == "function") cb = opts, opts = {};
  this.ondata = cb;
  return opts;
}
__name(StrmOpt, "StrmOpt");
class Deflate {
  static {
    __name(this, "Deflate");
  }
  constructor(opts, cb) {
    if (typeof opts == "function") cb = opts, opts = {};
    this.ondata = cb;
    this.o = opts || {};
    this.s = { l: 0, i: 32768, w: 32768, z: 32768 };
    this.b = new u8(98304);
    if (this.o.dictionary) {
      const dict = this.o.dictionary.subarray(-32768);
      this.b.set(dict, 32768 - dict.length);
      this.s.i = 32768 - dict.length;
    }
  }
  b;
  s;
  o;
  /**
   * The handler to call whenever data is available
   */
  ondata;
  p(c, f) {
    this.ondata(dopt(c, this.o, 0, 0, this.s), f);
  }
  /**
   * Pushes a chunk to be deflated
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    if (!this.ondata) err(5);
    if (this.s.l) err(4);
    const endLen = chunk.length + this.s.z;
    if (endLen > this.b.length) {
      if (endLen > 2 * this.b.length - 32768) {
        const newBuf = new u8(endLen & -32768);
        newBuf.set(this.b.subarray(0, this.s.z));
        this.b = newBuf;
      }
      const split = this.b.length - this.s.z;
      this.b.set(chunk.subarray(0, split), this.s.z);
      this.s.z = this.b.length;
      this.p(this.b, false);
      this.b.set(this.b.subarray(-32768));
      this.b.set(chunk.subarray(split), 32768);
      this.s.z = chunk.length - split + 32768;
      this.s.i = 32766, this.s.w = 32768;
    } else {
      this.b.set(chunk, this.s.z);
      this.s.z += chunk.length;
    }
    this.s.l = final & 1;
    if (this.s.z > this.s.w + 8191 || final) {
      this.p(this.b, final || false);
      this.s.w = this.s.i, this.s.i -= 2;
    }
  }
  /**
   * Flushes buffered uncompressed data. Useful to immediately retrieve the
   * deflated output for small inputs.
   */
  flush() {
    if (!this.ondata) err(5);
    if (this.s.l) err(4);
    this.p(this.b, false);
    this.s.w = this.s.i, this.s.i -= 2;
  }
}
class AsyncDeflate {
  static {
    __name(this, "AsyncDeflate");
  }
  /**
   * The handler to call whenever data is available
   */
  ondata;
  /**
   * The handler to call whenever buffered source data is processed (i.e. `queuedSize` updates)
   */
  ondrain;
  /**
   * The number of uncompressed bytes buffered in the stream
   */
  queuedSize;
  constructor(opts, cb) {
    astrmify([
      bDflt,
      () => [astrm, Deflate]
    ], this, StrmOpt.call(this, opts, cb), (ev) => {
      const strm = new Deflate(ev.data);
      onmessage = astrm(strm);
    }, 6, 1);
  }
  /**
   * A method to terminate the stream's internal worker. Subsequent calls to
   * push() will silently fail.
   */
  terminate;
}
function deflate(data, opts, cb) {
  if (!cb) cb = opts, opts = {};
  if (typeof cb != "function") err(7);
  return cbify(data, opts, [
    bDflt
  ], (ev) => pbf(deflateSync(ev.data[0], ev.data[1])), 0, cb);
}
__name(deflate, "deflate");
function deflateSync(data, opts) {
  return dopt(data, opts || {}, 0, 0);
}
__name(deflateSync, "deflateSync");
class Inflate {
  static {
    __name(this, "Inflate");
  }
  s;
  o;
  p;
  d;
  /**
   * The handler to call whenever data is available
   */
  ondata;
  constructor(opts, cb) {
    if (typeof opts == "function") cb = opts, opts = {};
    this.ondata = cb;
    const dict = opts && opts.dictionary && opts.dictionary.subarray(-32768);
    this.s = { i: 0, b: dict ? dict.length : 0 };
    this.o = new u8(32768);
    this.p = new u8(0);
    if (dict) this.o.set(dict);
  }
  e(c) {
    if (!this.ondata) err(5);
    if (this.d) err(4);
    if (!this.p.length) this.p = c;
    else if (c.length) {
      const n = new u8(this.p.length + c.length);
      n.set(this.p), n.set(c, this.p.length), this.p = n;
    }
  }
  c(final) {
    this.s.i = +(this.d = final || false);
    const bts = this.s.b;
    const dt = inflt(this.p, this.s, this.o);
    this.ondata(slc(dt, bts, this.s.b), this.d);
    this.o = slc(dt, this.s.b - 32768), this.s.b = this.o.length;
    this.p = slc(this.p, this.s.p / 8 | 0), this.s.p &= 7;
  }
  /**
   * Pushes a chunk to be inflated
   * @param chunk The chunk to push
   * @param final Whether this is the final chunk
   */
  push(chunk, final) {
    this.e(chunk), this.c(final);
  }
}
class AsyncInflate {
  static {
    __name(this, "AsyncInflate");
  }
  /**
   * The handler to call whenever data is available
   */
  ondata;
  /**
   * The handler to call whenever buffered source data is processed (i.e. `queuedSize` updates)
   */
  ondrain;
  /**
   * The number of compressed bytes buffered in the stream
   */
  queuedSize;
  constructor(opts, cb) {
    astrmify([
      bInflt,
      () => [astrm, Inflate]
    ], this, StrmOpt.call(this, opts, cb), (ev) => {
      const strm = new Inflate(ev.data);
      onmessage = astrm(strm);
    }, 7, 0);
  }
  /**
   * A method to terminate the stream's internal worker. Subsequent calls to
   * push() will silently fail.
   */
  terminate;
}
function inflate(data, opts, cb) {
  if (!cb) cb = opts, opts = {};
  if (typeof cb != "function") err(7);
  return cbify(data, opts, [
    bInflt
  ], (ev) => pbf(inflateSync(ev.data[0], gopt(ev.data[1]))), 1, cb);
}
__name(inflate, "inflate");
function inflateSync(data, opts) {
  return inflt(data, { i: 2 }, opts && opts.out, opts && opts.dictionary);
}
__name(inflateSync, "inflateSync");
class Gzip {
  static {
    __name(this, "Gzip");
  }
  c = crc();
  l = 0;
  v = 1;
  o;
  s;
  /**
   * The handler to call whenever data is available
   */
  ondata;
  constructor(opts, cb) {
    Deflate.call(this, opts, cb);
  }
  /**
   * Pushes a chunk to be GZIPped
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    this.c.p(chunk);
    this.l += chunk.length;
    Deflate.prototype.push.call(this, chunk, final);
  }
  p(c, f) {
    const raw = dopt(c, this.o, this.v && gzhl(this.o), f && 8, this.s);
    if (this.v) gzh(raw, this.o), this.v = 0;
    if (f) wbytes(raw, raw.length - 8, this.c.d()), wbytes(raw, raw.length - 4, this.l);
    this.ondata(raw, f);
  }
  /**
   * Flushes buffered uncompressed data. Useful to immediately retrieve the
   * GZIPped output for small inputs.
   */
  flush() {
    Deflate.prototype.flush.call(this);
  }
}
class AsyncGzip {
  static {
    __name(this, "AsyncGzip");
  }
  /**
   * The handler to call whenever data is available
   */
  ondata;
  /**
   * The handler to call whenever buffered source data is processed (i.e. `queuedSize` updates)
   */
  ondrain;
  /**
   * The number of uncompressed bytes buffered in the stream
   */
  queuedSize;
  constructor(opts, cb) {
    astrmify([
      bDflt,
      gze,
      () => [astrm, Deflate, Gzip]
    ], this, StrmOpt.call(this, opts, cb), (ev) => {
      const strm = new Gzip(ev.data);
      onmessage = astrm(strm);
    }, 8, 1);
  }
  /**
   * A method to terminate the stream's internal worker. Subsequent calls to
   * push() will silently fail.
   */
  terminate;
}
function gzip(data, opts, cb) {
  if (!cb) cb = opts, opts = {};
  if (typeof cb != "function") err(7);
  return cbify(data, opts, [
    bDflt,
    gze,
    () => [gzipSync]
  ], (ev) => pbf(gzipSync(ev.data[0], ev.data[1])), 2, cb);
}
__name(gzip, "gzip");
function gzipSync(data, opts) {
  if (!opts) opts = {};
  const c = crc(), l = data.length;
  c.p(data);
  const d = dopt(data, opts, gzhl(opts), 8), s = d.length;
  return gzh(d, opts), wbytes(d, s - 8, c.d()), wbytes(d, s - 4, l), d;
}
__name(gzipSync, "gzipSync");
class Gunzip {
  static {
    __name(this, "Gunzip");
  }
  v = 1;
  r = 0;
  o;
  p;
  s;
  /**
   * The handler to call whenever data is available
   */
  ondata;
  /**
   * The handler to call whenever a new GZIP member is found
   */
  onmember;
  constructor(opts, cb) {
    Inflate.call(this, opts, cb);
  }
  /**
   * Pushes a chunk to be GUNZIPped
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    Inflate.prototype.e.call(this, chunk);
    this.r += chunk.length;
    if (this.v) {
      const p = this.p.subarray(this.v - 1);
      const s = p.length > 3 ? gzs(p) : 4;
      if (s > p.length) {
        if (!final) return;
      } else if (this.v > 1 && this.onmember) {
        this.onmember(this.r - p.length);
      }
      this.p = p.subarray(s), this.v = 0;
    }
    Inflate.prototype.c.call(this, 0);
    if (this.s.f && !this.s.l) {
      this.v = shft(this.s.p) + 9;
      this.s = { i: 0 };
      this.o = new u8(0);
      this.push(new u8(0), final);
    } else if (final) {
      Inflate.prototype.c.call(this, final);
    }
  }
}
class AsyncGunzip {
  static {
    __name(this, "AsyncGunzip");
  }
  /**
   * The handler to call whenever data is available
   */
  ondata;
  /**
   * The handler to call whenever buffered source data is processed (i.e. `queuedSize` updates)
   */
  ondrain;
  /**
   * The number of compressed bytes buffered in the stream
   */
  queuedSize;
  /**
   * The handler to call whenever a new GZIP member is found
   */
  onmember;
  constructor(opts, cb) {
    astrmify([
      bInflt,
      guze,
      () => [astrm, Inflate, Gunzip]
    ], this, StrmOpt.call(this, opts, cb), (ev) => {
      const strm = new Gunzip(ev.data);
      strm.onmember = (offset) => postMessage(offset);
      onmessage = astrm(strm);
    }, 9, 0, (offset) => this.onmember && this.onmember(offset));
  }
  /**
   * A method to terminate the stream's internal worker. Subsequent calls to
   * push() will silently fail.
   */
  terminate;
}
function gunzip(data, opts, cb) {
  if (!cb) cb = opts, opts = {};
  if (typeof cb != "function") err(7);
  return cbify(data, opts, [
    bInflt,
    guze,
    () => [gunzipSync]
  ], (ev) => pbf(gunzipSync(ev.data[0], ev.data[1])), 3, cb);
}
__name(gunzip, "gunzip");
function gunzipSync(data, opts) {
  const st = gzs(data);
  if (st + 8 > data.length) err(6, "invalid gzip data");
  return inflt(data.subarray(st, -8), { i: 2 }, opts && opts.out || new u8(gzl(data)), opts && opts.dictionary);
}
__name(gunzipSync, "gunzipSync");
class Zlib {
  static {
    __name(this, "Zlib");
  }
  c = adler();
  v = 1;
  o;
  s;
  /**
   * The handler to call whenever data is available
   */
  ondata;
  constructor(opts, cb) {
    Deflate.call(this, opts, cb);
  }
  /**
   * Pushes a chunk to be zlibbed
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    this.c.p(chunk);
    Deflate.prototype.push.call(this, chunk, final);
  }
  p(c, f) {
    const raw = dopt(c, this.o, this.v && (this.o.dictionary ? 6 : 2), f && 4, this.s);
    if (this.v) zlh(raw, this.o), this.v = 0;
    if (f) wbytes(raw, raw.length - 4, this.c.d());
    this.ondata(raw, f);
  }
  /**
   * Flushes buffered uncompressed data. Useful to immediately retrieve the
   * zlibbed output for small inputs.
   */
  flush() {
    Deflate.prototype.flush.call(this);
  }
}
class AsyncZlib {
  static {
    __name(this, "AsyncZlib");
  }
  /**
   * The handler to call whenever data is available
   */
  ondata;
  /**
   * The handler to call whenever buffered source data is processed (i.e. `queuedSize` updates)
   */
  ondrain;
  /**
   * The number of uncompressed bytes buffered in the stream
   */
  queuedSize;
  constructor(opts, cb) {
    astrmify([
      bDflt,
      zle,
      () => [astrm, Deflate, Zlib]
    ], this, StrmOpt.call(this, opts, cb), (ev) => {
      const strm = new Zlib(ev.data);
      onmessage = astrm(strm);
    }, 10, 1);
  }
  /**
   * A method to terminate the stream's internal worker. Subsequent calls to
   * push() will silently fail.
   */
  terminate;
}
function zlib(data, opts, cb) {
  if (!cb) cb = opts, opts = {};
  if (typeof cb != "function") err(7);
  return cbify(data, opts, [
    bDflt,
    zle,
    () => [zlibSync]
  ], (ev) => pbf(zlibSync(ev.data[0], ev.data[1])), 4, cb);
}
__name(zlib, "zlib");
function zlibSync(data, opts) {
  if (!opts) opts = {};
  const a = adler();
  a.p(data);
  const d = dopt(data, opts, opts.dictionary ? 6 : 2, 4);
  return zlh(d, opts), wbytes(d, d.length - 4, a.d()), d;
}
__name(zlibSync, "zlibSync");
class Unzlib {
  static {
    __name(this, "Unzlib");
  }
  v;
  p;
  /**
   * The handler to call whenever data is available
   */
  ondata;
  constructor(opts, cb) {
    Inflate.call(this, opts, cb);
    this.v = opts && opts.dictionary ? 2 : 1;
  }
  /**
   * Pushes a chunk to be unzlibbed
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    Inflate.prototype.e.call(this, chunk);
    if (this.v) {
      if (this.p.length < 6 && !final) return;
      this.p = this.p.subarray(zls(this.p, this.v - 1)), this.v = 0;
    }
    if (final) {
      if (this.p.length < 4) err(6, "invalid zlib data");
      this.p = this.p.subarray(0, -4);
    }
    Inflate.prototype.c.call(this, final);
  }
}
class AsyncUnzlib {
  static {
    __name(this, "AsyncUnzlib");
  }
  /**
   * The handler to call whenever data is available
   */
  ondata;
  /**
   * The handler to call whenever buffered source data is processed (i.e. `queuedSize` updates)
   */
  ondrain;
  /**
   * The number of compressed bytes buffered in the stream
   */
  queuedSize;
  constructor(opts, cb) {
    astrmify([
      bInflt,
      zule,
      () => [astrm, Inflate, Unzlib]
    ], this, StrmOpt.call(this, opts, cb), (ev) => {
      const strm = new Unzlib(ev.data);
      onmessage = astrm(strm);
    }, 11, 0);
  }
  /**
   * A method to terminate the stream's internal worker. Subsequent calls to
   * push() will silently fail.
   */
  terminate;
}
function unzlib(data, opts, cb) {
  if (!cb) cb = opts, opts = {};
  if (typeof cb != "function") err(7);
  return cbify(data, opts, [
    bInflt,
    zule,
    () => [unzlibSync]
  ], (ev) => pbf(unzlibSync(ev.data[0], gopt(ev.data[1]))), 5, cb);
}
__name(unzlib, "unzlib");
function unzlibSync(data, opts) {
  return inflt(data.subarray(zls(data, opts && opts.dictionary), -4), { i: 2 }, opts && opts.out, opts && opts.dictionary);
}
__name(unzlibSync, "unzlibSync");
class Decompress {
  static {
    __name(this, "Decompress");
  }
  G;
  I;
  Z;
  o;
  s;
  p;
  /**
   * The handler to call whenever data is available
   */
  ondata;
  constructor(opts, cb) {
    this.o = StrmOpt.call(this, opts, cb) || {};
    this.G = Gunzip;
    this.I = Inflate;
    this.Z = Unzlib;
  }
  // init substream
  // overriden by AsyncDecompress
  i() {
    this.s.ondata = (dat, final) => {
      this.ondata(dat, final);
    };
  }
  /**
   * Pushes a chunk to be decompressed
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    if (!this.ondata) err(5);
    if (!this.s) {
      if (this.p && this.p.length) {
        const n = new u8(this.p.length + chunk.length);
        n.set(this.p), n.set(chunk, this.p.length);
      } else this.p = chunk;
      if (this.p.length > 2) {
        this.s = this.p[0] == 31 && this.p[1] == 139 && this.p[2] == 8 ? new this.G(this.o) : (this.p[0] & 15) != 8 || this.p[0] >> 4 > 7 || (this.p[0] << 8 | this.p[1]) % 31 ? new this.I(this.o) : new this.Z(this.o);
        this.i();
        this.s.push(this.p, final);
        this.p = null;
      }
    } else this.s.push(chunk, final);
  }
}
class AsyncDecompress {
  static {
    __name(this, "AsyncDecompress");
  }
  G;
  I;
  Z;
  /**
   * The handler to call whenever data is available
   */
  ondata;
  /**
   * The handler to call whenever buffered source data is processed (i.e. `queuedSize` updates)
   */
  ondrain;
  /**
   * The number of compressed bytes buffered in the stream
   */
  queuedSize;
  constructor(opts, cb) {
    Decompress.call(this, opts, cb);
    this.queuedSize = 0;
    this.G = AsyncGunzip;
    this.I = AsyncInflate;
    this.Z = AsyncUnzlib;
  }
  i() {
    this.s.ondata = (err2, dat, final) => {
      this.ondata(err2, dat, final);
    };
    this.s.ondrain = (size) => {
      this.queuedSize -= size;
      if (this.ondrain) this.ondrain(size);
    };
  }
  /**
   * Pushes a chunk to be decompressed
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    this.queuedSize += chunk.length;
    Decompress.prototype.push.call(this, chunk, final);
  }
}
function decompress(data, opts, cb) {
  if (!cb) cb = opts, opts = {};
  if (typeof cb != "function") err(7);
  return data[0] == 31 && data[1] == 139 && data[2] == 8 ? gunzip(data, opts, cb) : (data[0] & 15) != 8 || data[0] >> 4 > 7 || (data[0] << 8 | data[1]) % 31 ? inflate(data, opts, cb) : unzlib(data, opts, cb);
}
__name(decompress, "decompress");
function decompressSync(data, opts) {
  return data[0] == 31 && data[1] == 139 && data[2] == 8 ? gunzipSync(data, opts) : (data[0] & 15) != 8 || data[0] >> 4 > 7 || (data[0] << 8 | data[1]) % 31 ? inflateSync(data, opts) : unzlibSync(data, opts);
}
__name(decompressSync, "decompressSync");
const fltn = /* @__PURE__ */ __name((d, p, t, o) => {
  for (const k in d) {
    let val = d[k], n = p + k, op = o;
    if (Array.isArray(val)) op = mrg(o, val[1]), val = val[0];
    if (val instanceof u8) t[n] = [val, op];
    else {
      t[n += "/"] = [new u8(0), op];
      fltn(val, n, t, o);
    }
  }
}, "fltn");
const te = typeof TextEncoder != "undefined" && /* @__PURE__ */ new TextEncoder();
const td = typeof TextDecoder != "undefined" && /* @__PURE__ */ new TextDecoder();
let tds = 0;
try {
  td.decode(et, { stream: true });
  tds = 1;
} catch (e) {
}
const dutf8 = /* @__PURE__ */ __name((d) => {
  for (let r = "", i = 0; ; ) {
    let c = d[i++];
    const eb = (c > 127) + (c > 223) + (c > 239);
    if (i + eb > d.length) return { s: r, r: slc(d, i - 1) };
    if (!eb) r += String.fromCharCode(c);
    else if (eb == 3) {
      c = ((c & 15) << 18 | (d[i++] & 63) << 12 | (d[i++] & 63) << 6 | d[i++] & 63) - 65536, r += String.fromCharCode(55296 | c >> 10, 56320 | c & 1023);
    } else if (eb & 1) r += String.fromCharCode((c & 31) << 6 | d[i++] & 63);
    else r += String.fromCharCode((c & 15) << 12 | (d[i++] & 63) << 6 | d[i++] & 63);
  }
}, "dutf8");
class DecodeUTF8 {
  static {
    __name(this, "DecodeUTF8");
  }
  p;
  t;
  /**
   * Creates a UTF-8 decoding stream
   * @param cb The callback to call whenever data is decoded
   */
  constructor(cb) {
    this.ondata = cb;
    if (tds) this.t = new TextDecoder();
    else this.p = et;
  }
  /**
   * Pushes a chunk to be decoded from UTF-8 binary
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    if (!this.ondata) err(5);
    final = !!final;
    if (this.t) {
      this.ondata(this.t.decode(chunk, { stream: true }), final);
      if (final) {
        if (this.t.decode().length) err(8);
        this.t = null;
      }
      return;
    }
    if (!this.p) err(4);
    const dat = new u8(this.p.length + chunk.length);
    dat.set(this.p);
    dat.set(chunk, this.p.length);
    const { s, r } = dutf8(dat);
    if (final) {
      if (r.length) err(8);
      this.p = null;
    } else this.p = r;
    this.ondata(s, final);
  }
  /**
   * The handler to call whenever data is available
   */
  ondata;
}
class EncodeUTF8 {
  static {
    __name(this, "EncodeUTF8");
  }
  d;
  /**
   * Creates a UTF-8 decoding stream
   * @param cb The callback to call whenever data is encoded
   */
  constructor(cb) {
    this.ondata = cb;
  }
  /**
   * Pushes a chunk to be encoded to UTF-8
   * @param chunk The string data to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    if (!this.ondata) err(5);
    if (this.d) err(4);
    this.ondata(strToU8(chunk), this.d = final || false);
  }
  /**
   * The handler to call whenever data is available
   */
  ondata;
}
function strToU8(str, latin1) {
  if (latin1) {
    const ar2 = new u8(str.length);
    for (let i = 0; i < str.length; ++i) ar2[i] = str.charCodeAt(i);
    return ar2;
  }
  if (te) return te.encode(str);
  const l = str.length;
  let ar = new u8(str.length + (str.length >> 1));
  let ai = 0;
  const w = /* @__PURE__ */ __name((v) => {
    ar[ai++] = v;
  }, "w");
  for (let i = 0; i < l; ++i) {
    if (ai + 5 > ar.length) {
      const n = new u8(ai + 8 + (l - i << 1));
      n.set(ar);
      ar = n;
    }
    let c = str.charCodeAt(i);
    if (c < 128 || latin1) w(c);
    else if (c < 2048) w(192 | c >> 6), w(128 | c & 63);
    else if (c > 55295 && c < 57344)
      c = 65536 + (c & 1023 << 10) | str.charCodeAt(++i) & 1023, w(240 | c >> 18), w(128 | c >> 12 & 63), w(128 | c >> 6 & 63), w(128 | c & 63);
    else w(224 | c >> 12), w(128 | c >> 6 & 63), w(128 | c & 63);
  }
  return slc(ar, 0, ai);
}
__name(strToU8, "strToU8");
function strFromU8(dat, latin1) {
  if (latin1) {
    let r = "";
    for (let i = 0; i < dat.length; i += 16384)
      r += String.fromCharCode.apply(null, dat.subarray(i, i + 16384));
    return r;
  } else if (td) {
    return td.decode(dat);
  } else {
    const { s, r } = dutf8(dat);
    if (r.length) err(8);
    return s;
  }
}
__name(strFromU8, "strFromU8");
;
const dbf = /* @__PURE__ */ __name((l) => l == 1 ? 3 : l < 6 ? 2 : l == 9 ? 1 : 0, "dbf");
const slzh = /* @__PURE__ */ __name((d, b) => b + 30 + b2(d, b + 26) + b2(d, b + 28), "slzh");
const zh = /* @__PURE__ */ __name((d, b, z) => {
  const fnl = b2(d, b + 28), fn = strFromU8(d.subarray(b + 46, b + 46 + fnl), !(b2(d, b + 8) & 2048)), es = b + 46 + fnl, bs = b4(d, b + 20);
  const [sc, su, off] = z && bs == 4294967295 ? z64e(d, es) : [bs, b4(d, b + 24), b4(d, b + 42)];
  return [b2(d, b + 10), sc, su, fn, es + b2(d, b + 30) + b2(d, b + 32), off];
}, "zh");
const z64e = /* @__PURE__ */ __name((d, b) => {
  for (; b2(d, b) != 1; b += 4 + b2(d, b + 2)) ;
  return [b8(d, b + 12), b8(d, b + 4), b8(d, b + 20)];
}, "z64e");
const exfl = /* @__PURE__ */ __name((ex) => {
  let le = 0;
  if (ex) {
    for (const k in ex) {
      const l = ex[k].length;
      if (l > 65535) err(9);
      le += l + 4;
    }
  }
  return le;
}, "exfl");
const wzh = /* @__PURE__ */ __name((d, b, f, fn, u, c, ce, co) => {
  const fl2 = fn.length, ex = f.extra, col = co && co.length;
  let exl = exfl(ex);
  wbytes(d, b, ce != null ? 33639248 : 67324752), b += 4;
  if (ce != null) d[b++] = 20, d[b++] = f.os;
  d[b] = 20, b += 2;
  d[b++] = f.flag << 1 | (c < 0 && 8), d[b++] = u && 8;
  d[b++] = f.compression & 255, d[b++] = f.compression >> 8;
  const dt = new Date(f.mtime == null ? Date.now() : f.mtime), y = dt.getFullYear() - 1980;
  if (y < 0 || y > 119) err(10);
  wbytes(d, b, y << 25 | dt.getMonth() + 1 << 21 | dt.getDate() << 16 | dt.getHours() << 11 | dt.getMinutes() << 5 | dt.getSeconds() >> 1), b += 4;
  if (c != -1) {
    wbytes(d, b, f.crc);
    wbytes(d, b + 4, c < 0 ? -c - 2 : c);
    wbytes(d, b + 8, f.size);
  }
  wbytes(d, b + 12, fl2);
  wbytes(d, b + 14, exl), b += 16;
  if (ce != null) {
    wbytes(d, b, col);
    wbytes(d, b + 6, f.attrs);
    wbytes(d, b + 10, ce), b += 14;
  }
  d.set(fn, b);
  b += fl2;
  if (exl) {
    for (const k in ex) {
      const exf = ex[k], l = exf.length;
      wbytes(d, b, +k);
      wbytes(d, b + 2, l);
      d.set(exf, b + 4), b += 4 + l;
    }
  }
  if (col) d.set(co, b), b += col;
  return b;
}, "wzh");
const wzf = /* @__PURE__ */ __name((o, b, c, d, e) => {
  wbytes(o, b, 101010256);
  wbytes(o, b + 8, c);
  wbytes(o, b + 10, c);
  wbytes(o, b + 12, d);
  wbytes(o, b + 16, e);
}, "wzf");
class ZipPassThrough {
  static {
    __name(this, "ZipPassThrough");
  }
  filename;
  crc;
  size;
  compression;
  os;
  attrs;
  comment;
  extra;
  mtime;
  ondata;
  c;
  /**
   * Creates a pass-through stream that can be added to ZIP archives
   * @param filename The filename to associate with this data stream
   */
  constructor(filename) {
    this.filename = filename;
    this.c = crc();
    this.size = 0;
    this.compression = 0;
  }
  /**
   * Processes a chunk and pushes to the output stream. You can override this
   * method in a subclass for custom behavior, but by default this passes
   * the data through. You must call this.ondata(err, chunk, final) at some
   * point in this method.
   * @param chunk The chunk to process
   * @param final Whether this is the last chunk
   */
  process(chunk, final) {
    this.ondata(null, chunk, final);
  }
  /**
   * Pushes a chunk to be added. If you are subclassing this with a custom
   * compression algorithm, note that you must push data from the source
   * file only, pre-compression.
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    if (!this.ondata) err(5);
    this.c.p(chunk);
    this.size += chunk.length;
    if (final) this.crc = this.c.d();
    this.process(chunk, final || false);
  }
}
class ZipDeflate {
  static {
    __name(this, "ZipDeflate");
  }
  filename;
  crc;
  size;
  compression;
  flag;
  os;
  attrs;
  comment;
  extra;
  mtime;
  ondata;
  d;
  /**
   * Creates a DEFLATE stream that can be added to ZIP archives
   * @param filename The filename to associate with this data stream
   * @param opts The compression options
   */
  constructor(filename, opts) {
    if (!opts) opts = {};
    ZipPassThrough.call(this, filename);
    this.d = new Deflate(opts, (dat, final) => {
      this.ondata(null, dat, final);
    });
    this.compression = 8;
    this.flag = dbf(opts.level);
  }
  process(chunk, final) {
    try {
      this.d.push(chunk, final);
    } catch (e) {
      this.ondata(e, null, final);
    }
  }
  /**
   * Pushes a chunk to be deflated
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    ZipPassThrough.prototype.push.call(this, chunk, final);
  }
}
class AsyncZipDeflate {
  static {
    __name(this, "AsyncZipDeflate");
  }
  filename;
  crc;
  size;
  compression;
  flag;
  os;
  attrs;
  comment;
  extra;
  mtime;
  ondata;
  d;
  terminate;
  /**
   * Creates an asynchronous DEFLATE stream that can be added to ZIP archives
   * @param filename The filename to associate with this data stream
   * @param opts The compression options
   */
  constructor(filename, opts) {
    if (!opts) opts = {};
    ZipPassThrough.call(this, filename);
    this.d = new AsyncDeflate(opts, (err2, dat, final) => {
      this.ondata(err2, dat, final);
    });
    this.compression = 8;
    this.flag = dbf(opts.level);
    this.terminate = this.d.terminate;
  }
  process(chunk, final) {
    this.d.push(chunk, final);
  }
  /**
   * Pushes a chunk to be deflated
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    ZipPassThrough.prototype.push.call(this, chunk, final);
  }
}
class Zip {
  static {
    __name(this, "Zip");
  }
  u;
  d;
  /**
   * Creates an empty ZIP archive to which files can be added
   * @param cb The callback to call whenever data for the generated ZIP archive
   *           is available
   */
  constructor(cb) {
    this.ondata = cb;
    this.u = [];
    this.d = 1;
  }
  /**
   * Adds a file to the ZIP archive
   * @param file The file stream to add
   */
  add(file) {
    if (!this.ondata) err(5);
    if (this.d & 2) this.ondata(err(4 + (this.d & 1) * 8, 0, 1), null, false);
    else {
      const f = strToU8(file.filename), fl2 = f.length;
      const com = file.comment, o = com && strToU8(com);
      const u = fl2 != file.filename.length || o && com.length != o.length;
      const hl = fl2 + exfl(file.extra) + 30;
      if (fl2 > 65535) this.ondata(err(11, 0, 1), null, false);
      const header = new u8(hl);
      wzh(header, 0, file, f, u, -1);
      let chks = [header];
      const pAll = /* @__PURE__ */ __name(() => {
        for (const chk of chks) this.ondata(null, chk, false);
        chks = [];
      }, "pAll");
      let tr = this.d;
      this.d = 0;
      const ind = this.u.length;
      const uf = mrg(file, {
        f,
        u,
        o,
        t: /* @__PURE__ */ __name(() => {
          if (file.terminate) file.terminate();
        }, "t"),
        r: /* @__PURE__ */ __name(() => {
          pAll();
          if (tr) {
            const nxt = this.u[ind + 1];
            if (nxt) nxt.r();
            else this.d = 1;
          }
          tr = 1;
        }, "r")
      });
      let cl = 0;
      file.ondata = (err2, dat, final) => {
        if (err2) {
          this.ondata(err2, dat, final);
          this.terminate();
        } else {
          cl += dat.length;
          chks.push(dat);
          if (final) {
            const dd = new u8(16);
            wbytes(dd, 0, 134695760);
            wbytes(dd, 4, file.crc);
            wbytes(dd, 8, cl);
            wbytes(dd, 12, file.size);
            chks.push(dd);
            uf.c = cl, uf.b = hl + cl + 16, uf.crc = file.crc, uf.size = file.size;
            if (tr) uf.r();
            tr = 1;
          } else if (tr) pAll();
        }
      };
      this.u.push(uf);
    }
  }
  /**
   * Ends the process of adding files and prepares to emit the final chunks.
   * This *must* be called after adding all desired files for the resulting
   * ZIP file to work properly.
   */
  end() {
    if (this.d & 2) {
      this.ondata(err(4 + (this.d & 1) * 8, 0, 1), null, true);
      return;
    }
    if (this.d) this.e();
    else this.u.push({
      r: /* @__PURE__ */ __name(() => {
        if (!(this.d & 1)) return;
        this.u.splice(-1, 1);
        this.e();
      }, "r"),
      t: /* @__PURE__ */ __name(() => {
      }, "t")
    });
    this.d = 3;
  }
  e() {
    let bt = 0, l = 0, tl = 0;
    for (const f of this.u) tl += 46 + f.f.length + exfl(f.extra) + (f.o ? f.o.length : 0);
    const out = new u8(tl + 22);
    for (const f of this.u) {
      wzh(out, bt, f, f.f, f.u, -f.c - 2, l, f.o);
      bt += 46 + f.f.length + exfl(f.extra) + (f.o ? f.o.length : 0), l += f.b;
    }
    wzf(out, bt, this.u.length, tl, l);
    this.ondata(null, out, true);
    this.d = 2;
  }
  /**
   * A method to terminate any internal workers used by the stream. Subsequent
   * calls to add() will fail.
   */
  terminate() {
    for (const f of this.u) f.t();
    this.d = 2;
  }
  /**
   * The handler to call whenever data is available
   */
  ondata;
}
function zip(data, opts, cb) {
  if (!cb) cb = opts, opts = {};
  if (typeof cb != "function") err(7);
  const r = {};
  fltn(data, "", r, opts);
  const k = Object.keys(r);
  let lft = k.length, o = 0, tot = 0;
  const slft = lft, files = new Array(lft);
  const term = [];
  const tAll = /* @__PURE__ */ __name(() => {
    for (let i = 0; i < term.length; ++i) term[i]();
  }, "tAll");
  let cbd = /* @__PURE__ */ __name((a, b) => {
    mt(() => {
      cb(a, b);
    });
  }, "cbd");
  mt(() => {
    cbd = cb;
  });
  const cbf = /* @__PURE__ */ __name(() => {
    const out = new u8(tot + 22), oe = o, cdl = tot - o;
    tot = 0;
    for (let i = 0; i < slft; ++i) {
      const f = files[i];
      try {
        const l = f.c.length;
        wzh(out, tot, f, f.f, f.u, l);
        const badd = 30 + f.f.length + exfl(f.extra);
        const loc = tot + badd;
        out.set(f.c, loc);
        wzh(out, o, f, f.f, f.u, l, tot, f.m), o += 16 + badd + (f.m ? f.m.length : 0), tot = loc + l;
      } catch (e) {
        return cbd(e, null);
      }
    }
    wzf(out, o, files.length, cdl, oe);
    cbd(null, out);
  }, "cbf");
  if (!lft) cbf();
  for (let i = 0; i < slft; ++i) {
    const fn = k[i];
    const [file, p] = r[fn];
    const c = crc(), size = file.length;
    c.p(file);
    const f = strToU8(fn), s = f.length;
    const com = p.comment, m = com && strToU8(com), ms = m && m.length;
    const exl = exfl(p.extra);
    const compression = p.level == 0 ? 0 : 8;
    const cbl = /* @__PURE__ */ __name((e, d) => {
      if (e) {
        tAll();
        cbd(e, null);
      } else {
        const l = d.length;
        files[i] = mrg(p, {
          size,
          crc: c.d(),
          c: d,
          f,
          m,
          u: s != fn.length || m && com.length != ms,
          compression
        });
        o += 30 + s + exl + l;
        tot += 76 + 2 * (s + exl) + (ms || 0) + l;
        if (!--lft) cbf();
      }
    }, "cbl");
    if (s > 65535) cbl(err(11, 0, 1), null);
    if (!compression) cbl(null, file);
    else if (size < 16e4) {
      try {
        cbl(null, deflateSync(file, p));
      } catch (e) {
        cbl(e, null);
      }
    } else term.push(deflate(file, p, cbl));
  }
  return tAll;
}
__name(zip, "zip");
function zipSync(data, opts) {
  if (!opts) opts = {};
  const r = {};
  const files = [];
  fltn(data, "", r, opts);
  let o = 0;
  let tot = 0;
  for (const fn in r) {
    const [file, p] = r[fn];
    const compression = p.level == 0 ? 0 : 8;
    const f = strToU8(fn), s = f.length;
    const com = p.comment, m = com && strToU8(com), ms = m && m.length;
    const exl = exfl(p.extra);
    if (s > 65535) err(11);
    const d = compression ? deflateSync(file, p) : file, l = d.length;
    const c = crc();
    c.p(file);
    files.push(mrg(p, {
      size: file.length,
      crc: c.d(),
      c: d,
      f,
      m,
      u: s != fn.length || m && com.length != ms,
      o,
      compression
    }));
    o += 30 + s + exl + l;
    tot += 76 + 2 * (s + exl) + (ms || 0) + l;
  }
  const out = new u8(tot + 22), oe = o, cdl = tot - o;
  for (let i = 0; i < files.length; ++i) {
    const f = files[i];
    wzh(out, f.o, f, f.f, f.u, f.c.length);
    const badd = 30 + f.f.length + exfl(f.extra);
    out.set(f.c, f.o + badd);
    wzh(out, o, f, f.f, f.u, f.c.length, f.o, f.m), o += 16 + badd + (f.m ? f.m.length : 0);
  }
  wzf(out, o, files.length, cdl, oe);
  return out;
}
__name(zipSync, "zipSync");
class UnzipPassThrough {
  static {
    __name(this, "UnzipPassThrough");
  }
  static compression = 0;
  ondata;
  push(data, final) {
    this.ondata(null, data, final);
  }
}
class UnzipInflate {
  static {
    __name(this, "UnzipInflate");
  }
  static compression = 8;
  i;
  ondata;
  /**
   * Creates a DEFLATE decompression that can be used in ZIP archives
   */
  constructor() {
    this.i = new Inflate((dat, final) => {
      this.ondata(null, dat, final);
    });
  }
  push(data, final) {
    try {
      this.i.push(data, final);
    } catch (e) {
      this.ondata(e, null, final);
    }
  }
}
class AsyncUnzipInflate {
  static {
    __name(this, "AsyncUnzipInflate");
  }
  static compression = 8;
  i;
  ondata;
  terminate;
  /**
   * Creates a DEFLATE decompression that can be used in ZIP archives
   */
  constructor(_, sz) {
    if (sz < 32e4) {
      this.i = new Inflate((dat, final) => {
        this.ondata(null, dat, final);
      });
    } else {
      this.i = new AsyncInflate((err2, dat, final) => {
        this.ondata(err2, dat, final);
      });
      this.terminate = this.i.terminate;
    }
  }
  push(data, final) {
    if (this.i.terminate) data = slc(data, 0);
    this.i.push(data, final);
  }
}
class Unzip {
  static {
    __name(this, "Unzip");
  }
  d;
  c;
  p;
  k;
  o;
  /**
   * Creates a ZIP decompression stream
   * @param cb The callback to call whenever a file in the ZIP archive is found
   */
  constructor(cb) {
    this.onfile = cb;
    this.k = [];
    this.o = {
      0: UnzipPassThrough
    };
    this.p = et;
  }
  /**
   * Pushes a chunk to be unzipped
   * @param chunk The chunk to push
   * @param final Whether this is the last chunk
   */
  push(chunk, final) {
    if (!this.onfile) err(5);
    if (!this.p) err(4);
    if (this.c > 0) {
      const len = Math.min(this.c, chunk.length);
      const toAdd = chunk.subarray(0, len);
      this.c -= len;
      if (this.d) this.d.push(toAdd, !this.c);
      else this.k[0].push(toAdd);
      chunk = chunk.subarray(len);
      if (chunk.length) return this.push(chunk, final);
    } else {
      let f = 0, i = 0, is, buf;
      if (!this.p.length) buf = chunk;
      else if (!chunk.length) buf = this.p;
      else {
        buf = new u8(this.p.length + chunk.length);
        buf.set(this.p), buf.set(chunk, this.p.length);
      }
      const l = buf.length, oc = this.c, add = oc && this.d;
      for (; i < l - 4; ++i) {
        const sig = b4(buf, i);
        if (sig == 67324752) {
          f = 1, is = i;
          this.d = null;
          this.c = 0;
          const bf = b2(buf, i + 6), cmp = b2(buf, i + 8), u = bf & 2048, dd = bf & 8, fnl = b2(buf, i + 26), es = b2(buf, i + 28);
          if (l > i + 30 + fnl + es) {
            const chks = [];
            this.k.unshift(chks);
            f = 2;
            let sc = b4(buf, i + 18), su = b4(buf, i + 22);
            const fn = strFromU8(buf.subarray(i + 30, i += 30 + fnl), !u);
            if (sc == 4294967295) {
              [sc, su] = dd ? [-2] : z64e(buf, i);
            } else if (dd) sc = -1;
            i += es;
            this.c = sc;
            let d;
            const file = {
              name: fn,
              compression: cmp,
              start: /* @__PURE__ */ __name(() => {
                if (!file.ondata) err(5);
                if (!sc) file.ondata(null, et, true);
                else {
                  const ctr = this.o[cmp];
                  if (!ctr) file.ondata(err(14, "unknown compression type " + cmp, 1), null, false);
                  d = sc < 0 ? new ctr(fn) : new ctr(fn, sc, su);
                  d.ondata = (err2, dat, final2) => {
                    file.ondata(err2, dat, final2);
                  };
                  for (const dat of chks) d.push(dat, false);
                  if (this.k[0] == chks && this.c) this.d = d;
                  else d.push(et, true);
                }
              }, "start"),
              terminate: /* @__PURE__ */ __name(() => {
                if (d && d.terminate) d.terminate();
              }, "terminate")
            };
            if (sc >= 0) file.size = sc, file.originalSize = su;
            this.onfile(file);
          }
          break;
        } else if (oc) {
          if (sig == 134695760) {
            is = i += 12 + (oc == -2 && 8), f = 3, this.c = 0;
            break;
          } else if (sig == 33639248) {
            is = i -= 4, f = 3, this.c = 0;
            break;
          }
        }
      }
      this.p = et;
      if (oc < 0) {
        const dat = f ? buf.subarray(0, is - 12 - (oc == -2 && 8) - (b4(buf, is - 16) == 134695760 && 4)) : buf.subarray(0, i);
        if (add) add.push(dat, !!f);
        else this.k[+(f == 2)].push(dat);
      }
      if (f & 2) return this.push(buf.subarray(i), final);
      this.p = buf.subarray(i);
    }
    if (final) {
      if (this.c) err(13);
      this.p = null;
    }
  }
  /**
   * Registers a decoder with the stream, allowing for files compressed with
   * the compression type provided to be expanded correctly
   * @param decoder The decoder constructor
   */
  register(decoder) {
    this.o[decoder.compression] = decoder;
  }
  /**
   * The handler to call whenever a file is discovered
   */
  onfile;
}
const mt = typeof queueMicrotask == "function" ? queueMicrotask : typeof setTimeout == "function" ? setTimeout : (fn) => {
  fn();
};
function unzip(data, opts, cb) {
  if (!cb) cb = opts, opts = {};
  if (typeof cb != "function") err(7);
  const term = [];
  const tAll = /* @__PURE__ */ __name(() => {
    for (let i = 0; i < term.length; ++i) term[i]();
  }, "tAll");
  const files = {};
  let cbd = /* @__PURE__ */ __name((a, b) => {
    mt(() => {
      cb(a, b);
    });
  }, "cbd");
  mt(() => {
    cbd = cb;
  });
  let e = data.length - 22;
  for (; b4(data, e) != 101010256; --e) {
    if (!e || data.length - e > 65558) {
      cbd(err(13, 0, 1), null);
      return tAll;
    }
  }
  ;
  let lft = b2(data, e + 8);
  if (lft) {
    let c = lft;
    let o = b4(data, e + 16);
    let z = o == 4294967295 || c == 65535;
    if (z) {
      let ze = b4(data, e - 12);
      z = b4(data, ze) == 101075792;
      if (z) {
        c = lft = b4(data, ze + 32);
        o = b4(data, ze + 48);
      }
    }
    const fltr = opts && opts.filter;
    for (let i = 0; i < c; ++i) {
      const [c2, sc, su, fn, no, off] = zh(data, o, z), b = slzh(data, off);
      o = no;
      const cbl = /* @__PURE__ */ __name((e2, d) => {
        if (e2) {
          tAll();
          cbd(e2, null);
        } else {
          if (d) files[fn] = d;
          if (!--lft) cbd(null, files);
        }
      }, "cbl");
      if (!fltr || fltr({
        name: fn,
        size: sc,
        originalSize: su,
        compression: c2
      })) {
        if (!c2) cbl(null, slc(data, b, b + sc));
        else if (c2 == 8) {
          const infl = data.subarray(b, b + sc);
          if (su < 524288 || sc > 0.8 * su) {
            try {
              cbl(null, inflateSync(infl, { out: new u8(su) }));
            } catch (e2) {
              cbl(e2, null);
            }
          } else term.push(inflate(infl, { size: su }, cbl));
        } else cbl(err(14, "unknown compression type " + c2, 1), null);
      } else cbl(null, null);
    }
  } else cbd(null, {});
  return tAll;
}
__name(unzip, "unzip");
function unzipSync(data, opts) {
  const files = {};
  let e = data.length - 22;
  for (; b4(data, e) != 101010256; --e) {
    if (!e || data.length - e > 65558) err(13);
  }
  ;
  let c = b2(data, e + 8);
  if (!c) return {};
  let o = b4(data, e + 16);
  let z = o == 4294967295 || c == 65535;
  if (z) {
    let ze = b4(data, e - 12);
    z = b4(data, ze) == 101075792;
    if (z) {
      c = b4(data, ze + 32);
      o = b4(data, ze + 48);
    }
  }
  const fltr = opts && opts.filter;
  for (let i = 0; i < c; ++i) {
    const [c2, sc, su, fn, no, off] = zh(data, o, z), b = slzh(data, off);
    o = no;
    if (!fltr || fltr({
      name: fn,
      size: sc,
      originalSize: su,
      compression: c2
    })) {
      if (!c2) files[fn] = slc(data, b, b + sc);
      else if (c2 == 8) files[fn] = inflateSync(data.subarray(b, b + sc), { out: new u8(su) });
      else err(14, "unknown compression type " + c2);
    }
  }
  return files;
}
__name(unzipSync, "unzipSync");
export {
  AsyncGzip as AsyncCompress,
  AsyncDecompress,
  AsyncDeflate,
  AsyncGunzip,
  AsyncGzip,
  AsyncInflate,
  AsyncUnzipInflate,
  AsyncUnzlib,
  AsyncZipDeflate,
  AsyncZlib,
  Gzip as Compress,
  DecodeUTF8,
  Decompress,
  Deflate,
  EncodeUTF8,
  FlateErrorCode,
  Gunzip,
  Gzip,
  Inflate,
  Unzip,
  UnzipInflate,
  UnzipPassThrough,
  Unzlib,
  Zip,
  ZipDeflate,
  ZipPassThrough,
  Zlib,
  gzip as compress,
  gzipSync as compressSync,
  decompress,
  decompressSync,
  deflate,
  deflateSync,
  gunzip,
  gunzipSync,
  gzip,
  gzipSync,
  inflate,
  inflateSync,
  strFromU8,
  strToU8,
  unzip,
  unzipSync,
  unzlib,
  unzlibSync,
  zip,
  zipSync,
  zlib,
  zlibSync
};
//# sourceMappingURL=index.js.map
