/* Pure scoring functions: no DOM, no network. Same maths as clinical/score_patient.R. */
(function (root) {
  "use strict";

  // Convert the user's input (in the unit shown) to the model scale: x * scale, optionally ln().
  function toModel(biomarker, value) {
    var v = value * biomarker.scale;
    return biomarker.log ? Math.log(v) : v;
  }

  // Gompertz PhenoAge-type score: returns biological age in years.
  function phenoAge(model, x, age) {
    var xb = model.rate;
    Object.keys(model.betas).forEach(function (k) {
      xb += model.betas[k] * (k === "age" ? age : x[k]);
    });
    var M = 1 - Math.exp(model.m_n * Math.exp(xb) / model.m_d);
    return Math.log(model.ba_n * Math.log(1 - M)) / model.ba_d + model.ba_i;
  }

  // Piecewise-linear percentile from the published Z cut points (approximate, not an exact ECDF).
  function percentile(z, cuts) {
    var pts = [[-5, 0], [cuts.p10, 10], [cuts.p25, 25], [cuts.p50, 50], [cuts.p75, 75], [cuts.p90, 90], [5, 100]];
    if (z <= pts[0][0]) return 0;
    if (z >= pts[pts.length - 1][0]) return 100;
    for (var i = 1; i < pts.length; i++) {
      if (z <= pts[i][0]) {
        var a = pts[i - 1], b = pts[i];
        return a[1] + (b[1] - a[1]) * (z - a[0]) / (b[0] - a[0]);
      }
    }
    return NaN;
  }

  function band(z, cuts) {
    if (z < cuts.p10) return "P1";
    if (z < cuts.p25) return "P2";
    if (z < cuts.p75) return "P3";
    if (z < cuts.p90) return "P4";
    return "P5";
  }

  // entry: one registry record. inputs: {key: value in displayed unit}. sex: "F" | "M".
  function score(entry, inputs, age, sex) {
    var x = {};
    entry.clock.biomarkers.forEach(function (b) { x[b.key] = toModel(b, Number(inputs[b.key])); });
    var ba = phenoAge(entry.clock.model, x, age);
    var rc = entry.recentering;
    var expected = rc.intercept + rc.slope_age * age + rc.slope_male * (sex === "M" ? 1 : 0);
    var residual = ba - expected;
    var z = residual / rc.resid_sd;
    return { ba: ba, advancement: ba - age, expected: expected, residual: residual, z: z,
             percentile: percentile(z, rc.z_cutpoints), band: band(z, rc.z_cutpoints) };
  }

  var api = { toModel: toModel, phenoAge: phenoAge, percentile: percentile, band: band, score: score };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.BioAgeClock = api;
})(typeof self !== "undefined" ? self : this);
