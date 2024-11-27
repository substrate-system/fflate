var __defProp = Object.defineProperty;
var __name = (target, value) => __defProp(target, "name", { value, configurable: true });
const ch2 = {};
var worker_default = /* @__PURE__ */ __name((c, id, msg, transfer, cb) => {
  const w = new Worker(ch2[id] ||= URL.createObjectURL(
    new Blob([
      c + ';addEventListener("error",function(e){e=e.error;postMessage({$e$:[e.message,e.code,e.stack]})})'
    ], { type: "text/javascript" })
  ));
  w.onmessage = (e) => {
    const d = e.data, ed = d.$e$;
    if (ed) {
      const err = new Error(ed[0]);
      err["code"] = ed[1];
      err.stack = ed[2];
      cb(err, null);
    } else cb(null, d);
  };
  w.postMessage(msg, transfer);
  return w;
}, "default");
export {
  worker_default as default
};
//# sourceMappingURL=worker.js.map
