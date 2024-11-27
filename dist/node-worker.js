var __defProp = Object.defineProperty;
var __name = (target, value) => __defProp(target, "name", { value, configurable: true });
let Worker;
const workerAdd = ";var __w=require('worker_threads');__w.parentPort.on('message',function(m){onmessage({data:m})}),postMessage=function(m,t){__w.parentPort.postMessage(m,t)},close=process.exit;self=global";
try {
  Worker = require("worker_threads").Worker;
} catch (e) {
}
var node_worker_default = Worker ? (c, _, msg, transfer, cb) => {
  let done = false;
  const w = new Worker(c + workerAdd, { eval: true }).on("error", (e) => cb(e, null)).on("message", (m) => cb(null, m)).on("exit", (c2) => {
    if (c2 && !done) cb(new Error("exited with code " + c2), null);
  });
  w.postMessage(msg, transfer);
  w.terminate = () => {
    done = true;
    return Worker.prototype.terminate.call(w);
  };
  return w;
} : (_, __, ___, ____, cb) => {
  setImmediate(() => cb(new Error("async operations unsupported - update to Node 12+ (or Node 10-11 with the --experimental-worker CLI flag)"), null));
  const NOP = /* @__PURE__ */ __name(() => {
  }, "NOP");
  return {
    terminate: NOP,
    postMessage: NOP
  };
};
export {
  node_worker_default as default
};
//# sourceMappingURL=node-worker.js.map
