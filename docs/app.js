/* Renders the registry. All registry strings are written with textContent, never innerHTML. */
(function () {
  "use strict";
  var clocks = [];

  function el(tag, cls, text) {
    var n = document.createElement(tag);
    if (cls) n.className = cls;
    if (text !== undefined && text !== null) n.textContent = String(text);
    return n;
  }
  function fmt(x, d) { return typeof x === "number" ? x.toFixed(d) : "n/a"; }
  function int(x) { return typeof x === "number" ? x.toLocaleString("en-US") : "n/a"; }
  function badge(status) { return el("span", "badge " + status, status); }
  function kv(parent, k, v) { var r = el("div", "kv"); r.appendChild(el("span", null, k)); r.appendChild(el("span", null, v)); parent.appendChild(r); }

  function renderTable() {
    var tb = document.getElementById("rows");
    tb.textContent = "";
    clocks.forEach(function (c) {
      var tr = el("tr");
      tr.tabIndex = 0;
      var labCell = el("td"); labCell.appendChild(el("div", null, c.lab.name)); labCell.appendChild(el("div", "chip", c.lab.country));
      tr.appendChild(labCell);
      tr.appendChild(el("td", null, c.clock.name));
      var m = el("td"), chips = el("div", "chips");
      c.clock.biomarkers.forEach(function (b) { chips.appendChild(el("span", "chip", b.label)); });
      m.appendChild(chips); tr.appendChild(m);
      tr.appendChild(el("td", "num", int(c.performance.n_scored)));
      tr.appendChild(el("td", "num", fmt(c.performance.r_with_age, 3)));
      tr.appendChild(el("td", "num", fmt(c.performance.mean_ba_minus_ca_years, 2)));
      tr.appendChild(el("td", "num", fmt(c.performance.rmse_identity_years, 2)));
      var s = el("td"); s.appendChild(badge(c.status)); tr.appendChild(s);
      var open = function () { showDetail(c.id); };
      tr.addEventListener("click", open);
      tr.addEventListener("keydown", function (e) { if (e.key === "Enter") open(); });
      tb.appendChild(tr);
    });
  }

  var SVGNS = "http://www.w3.org/2000/svg";
  function renderMap() {
    var M = window.WORLD_MAP, svg = document.getElementById("arMap"), list = document.getElementById("mapList");
    if (!M || !svg) return;
    svg.setAttribute("viewBox", "0 0 " + M.w + " " + M.h);
    var land = document.createElementNS(SVGNS, "path");
    land.setAttribute("d", M.path); land.setAttribute("class", "land"); svg.appendChild(land);
    var pinned = clocks.filter(function (c) { return typeof c.lab.lat === "number" && typeof c.lab.lon === "number"; });
    pinned.forEach(function (c) {
      var x = (c.lab.lon - M.lon0) * M.k * M.s, y = (M.lat1 - c.lab.lat) * M.s;
      var dot = document.createElementNS(SVGNS, "circle");
      var pl = document.createElementNS(SVGNS, "circle");
      pl.setAttribute("cx", x.toFixed(1)); pl.setAttribute("cy", y.toFixed(1)); pl.setAttribute("r", 8); pl.setAttribute("class", "pulse");
      svg.appendChild(pl);
      dot.setAttribute("cx", x.toFixed(1)); dot.setAttribute("cy", y.toFixed(1)); dot.setAttribute("r", 8);
      dot.setAttribute("class", "pin" + (c.example ? " example" : "")); dot.setAttribute("tabindex", "0");
      dot.setAttribute("role", "button"); dot.setAttribute("aria-label", (c.lab.city || c.lab.name) + ": " + c.clock.name);
      var open = function () { showDetail(c.id); document.getElementById("detail").scrollIntoView({ behavior: "smooth", block: "nearest" }); };
      dot.addEventListener("click", open);
      dot.addEventListener("keydown", function (e) { if (e.key === "Enter") open(); });
      svg.appendChild(dot);
      var li = el("li"), b = el("button", null, (c.lab.city || c.lab.country) + " · " + c.clock.name);
      b.type = "button"; b.addEventListener("click", open); li.appendChild(b); list.appendChild(li);
    });
    if (!pinned.length) list.appendChild(el("li", null, "No entry has given a location yet."));
  }

  function renderStrip() {
    var box = document.getElementById("heroStats"); if (!box) return;
    var real = clocks.filter(function (c) { return !c.example; });
    var countries = {}; real.forEach(function (c) { countries[c.lab.country] = 1; });
    var scored = real.reduce(function (t, c) { return t + (c.performance.n_scored || 0); }, 0);
    [[real.length, real.length === 1 ? "registered clock" : "registered clocks"], [Object.keys(countries).length, Object.keys(countries).length === 1 ? "country" : "countries"], [int(scored), "people scored"]].forEach(function (x) {
      var d = el("div"); d.appendChild(el("b", null, x[0])); d.appendChild(el("span", null, x[1])); box.appendChild(d);
    });
  }

  function showDetail(id) {
    var c = clocks.filter(function (x) { return x.id === id; })[0];
    var box = document.getElementById("detail");
    if (!c) { box.classList.add("hidden"); return; }
    box.textContent = ""; box.classList.remove("hidden");
    var head = el("div"); head.appendChild(el("h2", null, c.clock.name + " · " + c.lab.name)); head.appendChild(badge(c.status));
    box.appendChild(head);
    box.appendChild(el("p", "sub", c.clock.description || ""));
    var g = el("div", "grid");
    var a = el("div", "card"); a.appendChild(el("h3", null, "Model"));
    kv(a, "Family", c.clock.family || "n/a"); kv(a, "Trained on", c.clock.trained_on || "n/a"); kv(a, "Evaluated on", c.clock.evaluated_on || "n/a");
    c.clock.biomarkers.forEach(function (b) { kv(a, b.label, b.unit + (b.log ? " (ln)" : "")); });
    var b2 = el("div", "card"); b2.appendChild(el("h3", null, "Local reference"));
    var r = c.recentering;
    kv(b2, "Method", r.method || "n/a"); kv(b2, "Reference n", int(r.reference_n));
    kv(b2, "Intercept", fmt(r.intercept, 3)); kv(b2, "Age slope", fmt(r.slope_age, 4)); kv(b2, "Male coefficient", fmt(r.slope_male, 3)); kv(b2, "Residual SD (y)", fmt(r.resid_sd, 3));
    ["p10", "p25", "p50", "p75", "p90"].forEach(function (k) { kv(b2, "Z " + k.toUpperCase(), fmt(r.z_cutpoints[k], 3)); });
    var p = el("div", "card"); p.appendChild(el("h3", null, "Performance"));
    kv(p, "n scored", int(c.performance.n_scored)); kv(p, "r with age", fmt(c.performance.r_with_age, 3));
    kv(p, "RMSE regression (y)", fmt(c.performance.rmse_regression_years, 2)); kv(p, "RMSE identity (y)", fmt(c.performance.rmse_identity_years, 2));
    kv(p, "Mean BA − CA (y)", fmt(c.performance.mean_ba_minus_ca_years, 2));
    var q = el("div", "card"); q.appendChild(el("h3", null, "Cohort"));
    var h = c.cohort || {};
    kv(q, "Eligible n", int(h.n_eligible)); kv(q, "Female (%)", fmt(h.pct_female, 1)); kv(q, "Mean age (y)", fmt(h.mean_age, 1)); kv(q, "License", c.license);
    [a, b2, p, q].forEach(function (x) { g.appendChild(x); });
    box.appendChild(g);
    var sel = document.getElementById("clockSel"); sel.value = c.id; buildFields();
    box.scrollIntoView({ behavior: "smooth", block: "nearest" });
  }

  function current() { var id = document.getElementById("clockSel").value; return clocks.filter(function (x) { return x.id === id; })[0]; }

  function buildFields() {
    var c = current(), f = document.getElementById("fields");
    f.textContent = ""; document.getElementById("out").classList.add("hidden"); document.getElementById("calcErr").textContent = "";
    function field(id, label, type, extra) {
      var w = el("div"); var l = el("label", null, label); l.htmlFor = id;
      var i = el("input"); i.id = id; i.type = type; i.inputMode = "decimal"; i.step = "any"; if (extra) Object.keys(extra).forEach(function (k) { i[k] = extra[k]; });
      w.appendChild(l); w.appendChild(i); f.appendChild(w); return i;
    }
    field("in_age", "Age (years)", "number", { min: 18, max: 120 });
    var sw = el("div"); var sl = el("label", null, "Sex"); sl.htmlFor = "in_sex"; var ss = el("select"); ss.id = "in_sex";
    [["F", "Female"], ["M", "Male"]].forEach(function (o) { var op = el("option", null, o[1]); op.value = o[0]; ss.appendChild(op); });
    sw.appendChild(sl); sw.appendChild(ss); f.appendChild(sw);
    c.clock.biomarkers.forEach(function (b) { field("in_" + b.key, b.label + " (" + b.unit + ")", "number", { min: 0 }); });
  }

  function calculate() {
    var c = current(), err = document.getElementById("calcErr"); err.textContent = "";
    var age = parseFloat(document.getElementById("in_age").value), sex = document.getElementById("in_sex").value, inputs = {};
    if (!(age >= 18 && age <= 120)) { err.textContent = "Enter an age between 18 and 120."; return; }
    for (var i = 0; i < c.clock.biomarkers.length; i++) {
      var b = c.clock.biomarkers[i], v = parseFloat(document.getElementById("in_" + b.key).value);
      if (!(v > 0)) { err.textContent = "Enter a positive value for " + b.label + "."; return; }
      inputs[b.key] = v;
    }
    var r = BioAgeClock.score(c, inputs, age, sex);
    if (!isFinite(r.ba)) { err.textContent = "These values fall outside what the model can compute."; return; }
    var s = document.getElementById("stats"); s.textContent = "";
    [["Biological age", fmt(r.ba, 1) + " y"], ["BA − CA", (r.advancement > 0 ? "+" : "") + fmt(r.advancement, 1) + " y"],
     ["Expected BA, same age and sex", fmt(r.expected, 1) + " y"], ["Local Z score", (r.z > 0 ? "+" : "") + fmt(r.z, 2)],
     ["Approx. percentile", fmt(r.percentile, 0)], ["Group", r.band]].forEach(function (x) {
      var d = el("div", "stat"); d.appendChild(el("div", "v", x[1])); d.appendChild(el("div", "l", x[0])); s.appendChild(d);
    });
    var bar = document.getElementById("bar"); bar.textContent = "";
    [["b1", 10, "P1"], ["b2", 15, "P2"], ["b3", 50, "P3"], ["b4", 15, "P4"], ["b5", 10, "P5"]].forEach(function (x) {
      var d = el("div", x[0], x[2]); d.style.width = x[1] + "%"; bar.appendChild(d);
    });
    var m = el("div", "marker"); m.style.left = "calc(" + Math.max(0, Math.min(100, r.percentile)) + "% - 1px)"; bar.appendChild(m);
    document.getElementById("out").classList.remove("hidden");
  }

  fetch("registry.json").then(function (r) { return r.json(); }).then(function (data) {
    clocks = data.clocks.slice().sort(function (a, b) { return (a.example === true) - (b.example === true); });
    renderTable();
    renderMap();
    renderStrip();
    var sel = document.getElementById("clockSel");
    clocks.forEach(function (c) { var o = el("option", null, c.clock.name + " · " + c.lab.name); o.value = c.id; sel.appendChild(o); });
    sel.addEventListener("change", buildFields); buildFields();
    document.getElementById("go").addEventListener("click", calculate);
    var m = /#clock=([a-z0-9-]+)/.exec(location.hash); if (m) showDetail(m[1]);
  }).catch(function () {
    var tb = document.getElementById("rows"); var tr = el("tr"); var td = el("td", null, "The registry could not be loaded. If you opened this file directly, serve the folder over HTTP (for example: python -m http.server)."); td.colSpan = 8; tr.appendChild(td); tb.appendChild(tr);
  });
})();
