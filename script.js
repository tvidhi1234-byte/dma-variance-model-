
const $ = id => document.getElementById(id);

const names = ["Method A", "Method B", "Method C"];
const colors = ["#8574df", "#e69bb8", "#65bdb4"];
const alpha = 0.05;

function parseGroup(id) {
  const values = $(id).value
    .trim()
    .split(/[,\s;]+/)
    .filter(Boolean)
    .map(Number);

  if (values.length < 2 ||
      values.some(v => !Number.isFinite(v))) {
    throw new Error(
      "Enter at least two valid numbers in every group."
    );
  }

  return values;
}

function mean(a) {
  return a.reduce((s, x) => s + x, 0) / a.length;
}

function variance(a) {
  const m = mean(a);
  return a.reduce((s, x) => s + (x - m) ** 2, 0)
    / (a.length - 1);
}

function median(a) {
  const s = [...a].sort((x, y) => x - y);
  const m = Math.floor(s.length / 2);

  return s.length % 2
    ? s[m]
    : (s[m - 1] + s[m]) / 2;
}

// Log-gamma approximation for statistical calculations.
function logGamma(z) {
  const p = [
    676.5203681218851, -1259.1392167224028,
    771.3234287776531, -176.6150291621406,
    12.507343278686905, -0.13857109526572012,
    9.984369578019572e-6, 1.5056327351493116e-7
  ];

  if (z < 0.5) {
    return Math.log(Math.PI)
      - Math.log(Math.sin(Math.PI * z))
      - logGamma(1 - z);
  }

  z -= 1;
  let x = 0.99999999999980993;

  p.forEach((v, i) => {
    x += v / (z + i + 1);
  });

  const t = z + p.length - 0.5;

  return 0.5 * Math.log(2 * Math.PI)
    + (z + 0.5) * Math.log(t)
    - t + Math.log(x);
}

function betaFraction(a, b, x) {
  const max = 200;
  const eps = 3e-14;
  const tiny = 1e-300;

  const qab = a + b;
  const qap = a + 1;
  const qam = a - 1;

  let c = 1;
  let d = 1 - qab * x / qap;

  if (Math.abs(d) < tiny) d = tiny;
  d = 1 / d;

  let h = d;

  for (let m = 1; m <= max; m++) {
    const m2 = 2 * m;

    let aa = m * (b - m) * x
      / ((qam + m2) * (a + m2));

    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    h *= d * c;

    aa = -(a + m) * (qab + m) * x
      / ((a + m2) * (qap + m2));

    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;

    const delta = d * c;
    h *= delta;

    if (Math.abs(delta - 1) < eps) break;
  }

  return h;
}

function beta(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;

  const bt = Math.exp(
    logGamma(a + b) - logGamma(a) - logGamma(b)
    + a * Math.log(x) + b * Math.log1p(-x)
  );

  if (x < (a + 1) / (a + b + 2)) {
    return bt * betaFraction(a, b, x) / a;
  }

  return 1 - bt * betaFraction(b, a, 1 - x) / b;
}

function fPValue(f, df1, df2) {
  if (!Number.isFinite(f)) return 0;
  if (f <= 0) return 1;

  return beta(
    df2 / (df2 + df1 * f),
    df2 / 2,
    df1 / 2
  );
}

// One-way ANOVA calculation.
function anova(groups) {
  const all = groups.flat();
  const grandMean = mean(all);
  const k = groups.length;
  const n = all.length;

  const ssBetween = groups.reduce((sum, g) =>
    sum + g.length * (mean(g) - grandMean) ** 2, 0);

  const ssWithin = groups.reduce((sum, g) =>
    sum + g.reduce((s, x) =>
      s + (x - mean(g)) ** 2, 0), 0);

  const df1 = k - 1;
  const df2 = n - k;

  const msBetween = ssBetween / df1;
  const msWithin = ssWithin / df2;

  const f = msWithin === 0
    ? (msBetween === 0 ? 0 : Infinity)
    : msBetween / msWithin;

  return {
    f,
    p: fPValue(f, df1, df2),
    df1,
    df2
  };
}

// Levene's test using absolute deviations from group medians.
function levene(groups) {
  const deviations = groups.map(g => {
    const med = median(g);
    return g.map(x => Math.abs(x - med));
  });

  return anova(deviations);
}

function fmt(n) {
  if (!Number.isFinite(n)) return "∞";
  if (n > 0 && n < 0.001) return n.toExponential(2);
  return n.toFixed(3);
}

function analyze() {
  $("error").textContent = "";

  try {
    const groups = [
      parseGroup("groupA"),
      parseGroup("groupB"),
      parseGroup("groupC")
    ];

    const a = anova(groups);
    const l = levene(groups);

    $("anovaP").textContent = fmt(a.p);
    $("leveneP").textContent = fmt(l.p);

    const aSig = a.p < alpha;
    const lSig = l.p < alpha;

    $("anovaStatus").textContent = aSig
      ? "Significant difference"
      : "No significant difference";

    $("leveneStatus").textContent = lSig
      ? "Variability differs"
      : "Similar variability";

    $("anovaExplanation").textContent = aSig
      ? `ANOVA: F(${a.df1}, ${a.df2}) = ${fmt(a.f)}. `
        + "The p-value is below 0.05. At least one group mean differs."
      : `ANOVA: F(${a.df1}, ${a.df2}) = ${fmt(a.f)}. `
        + "There is insufficient evidence of a difference in group means.";

    $("leveneExplanation").textContent = lSig
      ? `Levene's test: F(${l.df1}, ${l.df2}) = ${fmt(l.f)}. `
        + "The equal-variance assumption may not hold."
      : `Levene's test: F(${l.df1}, ${l.df2}) = ${fmt(l.f)}. `
        + "There is insufficient evidence that group variances differ.";

    // Draw average-mark bars.
    const means = groups.map(mean);
    const maxMean = Math.max(...means.map(Math.abs), 1);

    $("chart").innerHTML = groups.map((g, i) => `
      <div class="bar-group">
        <div class="bar-value">${fmt(means[i])}</div>
        <div class="bar"
          style="height:${Math.max(5, Math.abs(means[i]) / maxMean * 125)}px;
          background:${colors[i]}"></div>
        <div class="bar-label">${names[i]}</div>
      </div>
    `).join("");

    // Display group statistics.
    $("summaryBody").innerHTML = groups.map((g, i) => `
      <tr>
        <td>${names[i]}</td>
        <td>${g.length}</td>
        <td>${fmt(mean(g))}</td>
        <td>${fmt(variance(g))}</td>
      </tr>
    `).join("");

  } catch (error) {
    $("error").textContent = error.message;
  }
}

$("analyzeBtn").addEventListener("click", analyze);
analyze();
