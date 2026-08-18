"use strict";

const el = (id) => document.getElementById(id);

el("tips-btn").addEventListener("click", () => {
  el("tips-panel").classList.toggle("hidden");
});

// Session-wide: once enabled, every study submitted from here on uses
// "pbs" as its executor instead of "local", until toggled back off.
let pbsEnabled = false;
let pbsTemplateLoaded = false;

el("pbs-toggle-btn").addEventListener("click", async () => {
  pbsEnabled = !pbsEnabled;
  el("pbs-toggle-btn").textContent = pbsEnabled ? "PBS submission: ON" : "Enable PBS submission";
  el("pbs-toggle-btn").classList.toggle("active", pbsEnabled);
  el("pbs-warning").classList.toggle("hidden", !pbsEnabled);

  // Fetch the default script once, the first time PBS mode is turned
  // on -- after that the box keeps whatever the user has edited, even
  // across toggling off and back on, so an in-progress edit is never
  // lost.
  if (pbsEnabled && !pbsTemplateLoaded) {
    const textarea = el("pbs-script-textarea");
    textarea.value = "Loading...";
    try {
      const res = await fetch("/api/pbs_template");
      const data = await res.json();
      textarea.value = data.template;
      pbsTemplateLoaded = true;
    } catch (e) {
      textarea.value = `Couldn't load the default template: ${e.message}`;
    }
  }
});

// The script starts read-only (view mode) so it can't be changed by
// accident -- "Edit" unlocks it, "Save" locks it back. The textarea's
// value is what's actually submitted either way; Save is a UI
// affordance confirming the edit is done, not a separate persistence
// step.
el("pbs-edit-btn").addEventListener("click", () => {
  const textarea = el("pbs-script-textarea");
  textarea.readOnly = false;
  textarea.focus();
  el("pbs-edit-btn").classList.add("hidden");
  el("pbs-save-btn").classList.remove("hidden");
});

el("pbs-save-btn").addEventListener("click", () => {
  el("pbs-script-textarea").readOnly = true;
  el("pbs-save-btn").classList.add("hidden");
  el("pbs-edit-btn").classList.remove("hidden");
});

let currentUpload = null; // { token, elements } for whichever geometry source is active

function geometrySource() {
  return document.querySelector('input[name="geometry_source"]:checked').value;
}

function updateGeometrySourceVisibility() {
  const source = geometrySource();
  el("field-gaussian-log").classList.toggle("hidden", source !== "gaussian_log");
  el("field-gamess-inp").classList.toggle("hidden", source !== "gamess_inp");
  el("field-guess").classList.toggle("hidden", source !== "guess");
}

document.querySelectorAll('input[name="geometry_source"]').forEach((r) =>
  r.addEventListener("change", () => {
    currentUpload = null;
    updateGeometrySourceVisibility();
  })
);
updateGeometrySourceVisibility();

// ---------------------------------------------------------------- uploads

async function postJSON(url, body) {
  const res = await fetch(url, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!res.ok) {
    const errBody = await res.json().catch(() => ({}));
    return { ok: false, message: errBody.detail || `HTTP ${res.status}` };
  }
  return { ok: true, body: await res.json() };
}

function applyUploadResult(infoId, data, errorPrefix) {
  const info = el(infoId);
  if (data.ok) {
    currentUpload = data.body;
    info.textContent = data.body.elements ? `OK: elements ${data.body.elements.join(", ")}` : "OK";
    info.className = "file-info ok";
  } else {
    currentUpload = null;
    info.textContent = `${errorPrefix}: ${data.message}`;
    info.className = "file-info err";
  }
}

async function handleFileUpload(inputId, infoId, kind) {
  const input = el(inputId);
  const info = el(infoId);
  const file = input.files[0];
  if (!file) {
    currentUpload = null;
    info.textContent = "";
    info.className = "file-info";
    return;
  }
  info.textContent = "Reading...";
  info.className = "file-info";

  const form = new FormData();
  form.append("file", file);
  let data;
  try {
    const res = await fetch(`/api/upload?kind=${kind}`, { method: "POST", body: form });
    if (!res.ok) {
      const errBody = await res.json().catch(() => ({}));
      data = { ok: false, message: errBody.detail || `HTTP ${res.status}` };
    } else {
      data = { ok: true, body: await res.json() };
    }
  } catch (e) {
    data = { ok: false, message: e.message };
  }
  applyUploadResult(infoId, data, "Couldn't read this file");
}

el("gaussian-log-file").addEventListener("change", () => handleFileUpload("gaussian-log-file", "gaussian-log-info", "gaussian_log"));
el("gamess-inp-file").addEventListener("change", () => handleFileUpload("gamess-inp-file", "gamess-inp-info", "gamess_inp"));
el("guess-file").addEventListener("change", () => handleFileUpload("guess-file", "guess-info", "guess"));

el("guess-smiles-load").addEventListener("click", async () => {
  const smiles = el("guess-smiles").value.trim();
  const info = el("guess-info");
  if (!smiles) {
    info.textContent = "Enter a SMILES string first.";
    info.className = "file-info err";
    return;
  }
  const button = el("guess-smiles-load");
  button.disabled = true;
  info.textContent = "Embedding...";
  info.className = "file-info";
  el("guess-file").value = ""; // a SMILES load supersedes any file already chosen

  const data = await postJSON("/api/analyze_smiles", { smiles, label: "guess" });
  applyUploadResult("guess-info", data, `Couldn't embed ${smiles}`);
  button.disabled = false;
});

el("guess-smiles").addEventListener("keydown", (evt) => {
  if (evt.key === "Enter") {
    evt.preventDefault();
    el("guess-smiles-load").click();
  }
});

// -------------------------------------------------------------- basis set

(async () => {
  const data = await fetch("/api/gbasis_choices").then((r) => r.json());
  const select = el("gbasis");
  Object.entries(data.choices).forEach(([label, line]) => {
    const opt = document.createElement("option");
    opt.value = line;
    opt.textContent = label;
    select.appendChild(opt);
  });
  const custom = document.createElement("option");
  custom.value = "__custom__";
  custom.textContent = "Custom...";
  select.appendChild(custom);
})();

el("gbasis").addEventListener("change", () => {
  el("gbasis-custom").classList.toggle("hidden", el("gbasis").value !== "__custom__");
});

// ------------------------------------ guess-geometry Gaussian level

(async () => {
  const data = await fetch("/api/gaussian_levels").then((r) => r.json());
  const methodSel = el("gaussian-method");
  data.methods.forEach((m) => {
    const opt = document.createElement("option");
    opt.value = m;
    opt.textContent = m;
    if (m === "B3LYP") opt.selected = true;
    methodSel.appendChild(opt);
  });
  const basisSel = el("gaussian-basis");
  data.basis_sets.forEach((b) => {
    const opt = document.createElement("option");
    opt.value = b;
    opt.textContent = b;
    if (b === "6-31G(d)") opt.selected = true;
    basisSel.appendChild(opt);
  });
})();

// --------------------------------------------------------- GAMESS settings

(async () => {
  const data = await fetch("/api/gamess_defaults").then((r) => r.json());
  el("rungms-path").value = data.rungms_path;
  el("scratch-dir").value = data.scratch_dir;
})();

// ------------------------------------------------------------------- CIS

function updateCisCasscfVisibility() {
  const cisOn = el("run-cis").checked;
  el("nstate-field").classList.toggle("hidden", !cisOn);
  el("run-casscf").disabled = !cisOn;
  if (!cisOn) {
    el("run-casscf").checked = false;
    el("casscf-fields").classList.add("hidden");
  }
}

el("run-cis").addEventListener("change", updateCisCasscfVisibility);
updateCisCasscfVisibility();

el("run-casscf").addEventListener("change", () => {
  el("casscf-fields").classList.toggle("hidden", !el("run-casscf").checked);
});

// ----------------------------------------------------------------- submit

function showError(message) {
  el("form-error").textContent = message;
}

el("job-form").addEventListener("submit", async (evt) => {
  evt.preventDefault();
  showError("");

  if (!currentUpload) {
    showError("Provide a geometry first (upload a file or, for a guess geometry, load a SMILES string).");
    return;
  }

  const gbasisSel = el("gbasis").value;
  const gbasisLine = gbasisSel === "__custom__" ? el("gbasis-custom").value.trim() : gbasisSel;

  const req = {
    job_name: el("job-name").value.trim() || "job",
    geometry_source: geometrySource(),
    upload_token: currentUpload.token,
    nprocs: parseInt(el("nprocs").value, 10) || 2,
    mem_gb: parseInt(el("mem-gb").value, 10) || 2,
    gaussian_method: el("gaussian-method").value,
    gaussian_basis: el("gaussian-basis").value,
    charge: parseInt(el("charge").value, 10) || 0,
    mult: parseInt(el("mult").value, 10) || 1,
    gbasis_line: gbasisLine,
    use_soscf: document.querySelector('input[name="soscf"]:checked').value === "soscf",
    run_cis: el("run-cis").checked,
    nstate: parseInt(el("nstate").value, 10) || 5,
    run_casscf: el("run-casscf").checked,
    casscf_threshold: parseFloat(el("casscf-threshold").value) || 0.20,
    rungms_path: el("rungms-path").value.trim(),
    scratch_dir: el("scratch-dir").value.trim(),
    ncpus: parseInt(el("ncpus").value, 10) || 1,
    mem_mwords: parseInt(el("mem-mwords").value, 10) || 1,
    executor: pbsEnabled ? "pbs" : "local",
    pbs_script_template: pbsEnabled ? el("pbs-script-textarea").value : null,
  };

  el("submit-btn").disabled = true;
  try {
    const res = await fetch("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(req),
    });
    if (!res.ok) {
      const body = await res.json().catch(() => ({}));
      throw new Error(body.detail || `HTTP ${res.status}`);
    }
    const { job_id } = await res.json();
    addStudyCard(job_id, req.job_name, parseInt(el("casscf-nstate").value, 10) || 3);
    el("job-form").classList.add("hidden");
  } catch (e) {
    showError(`Couldn't start the job: ${e.message}`);
  } finally {
    el("submit-btn").disabled = false;
  }
});

el("add-study-btn").addEventListener("click", () => {
  showError("");
  el("job-form").classList.remove("hidden");
  el("job-form").scrollIntoView({ behavior: "smooth", block: "start" });
});

// ---------------------------------------------------------------- studies

const TRASH_ICON = `<svg viewBox="0 0 16 16" width="14" height="14" fill="none" stroke="currentColor" stroke-width="1.4" stroke-linecap="round" stroke-linejoin="round"><path d="M2 4h12M5.5 4V2.5A1 1 0 0 1 6.5 1.5h3A1 1 0 0 1 10.5 2.5V4M6.5 7.5v4.5M9.5 7.5v4.5M3.5 4l.6 8.4A1 1 0 0 0 5.1 13.4h5.8a1 1 0 0 0 1-1.1L13 4"/></svg>`;

function addStudyCard(jobId, jobName, defaultCasscfNstate) {
  const list = el("studies-list");
  const card = document.createElement("div");
  card.className = "study-card";
  card.dataset.jobId = jobId;
  card.innerHTML = `
    <div class="study-card-header">
      <span class="study-name">${jobName}</span>
      <div class="study-card-actions">
        <button type="button" class="stop-btn small-btn">Stop</button>
        <span class="badge running">Running</span>
        <button type="button" class="delete-btn" title="Remove this study entry" disabled>${TRASH_ICON}</button>
      </div>
    </div>
    <p class="stage-text"></p>
    <details class="study-log-details">
      <summary>Log</summary>
      <pre class="log-panel"></pre>
    </details>
    <div class="study-results"></div>
  `;
  list.prepend(card);
  el("studies-panel").classList.remove("hidden");

  const deleteBtn = card.querySelector(".delete-btn");
  deleteBtn.addEventListener("click", async () => {
    deleteBtn.disabled = true;
    try {
      const res = await fetch(`/api/jobs/${jobId}`, { method: "DELETE" });
      if (!res.ok) {
        const body = await res.json().catch(() => ({}));
        throw new Error(body.detail || `HTTP ${res.status}`);
      }
      card.remove();
      if (!list.children.length) {
        el("studies-panel").classList.add("hidden");
        el("job-form").classList.remove("hidden");
      }
    } catch (e) {
      deleteBtn.disabled = false;
      deleteBtn.title = `Couldn't remove: ${e.message}`;
    }
  });

  const stopBtn = card.querySelector(".stop-btn");
  stopBtn.addEventListener("click", async () => {
    stopBtn.disabled = true;
    stopBtn.textContent = "Stopping...";
    try {
      await fetch(`/api/jobs/${jobId}/cancel`, { method: "POST" });
    } catch (e) {
      stopBtn.disabled = false;
      stopBtn.textContent = "Stop";
    }
  });

  pollStudy(jobId, card, defaultCasscfNstate);
}

function pollStudy(jobId, card, defaultCasscfNstate) {
  const badge = card.querySelector(".badge");
  const stopBtn = card.querySelector(".stop-btn");
  const deleteBtn = card.querySelector(".delete-btn");
  const stageText = card.querySelector(".stage-text");
  const logPanel = card.querySelector(".log-panel");
  const resultsDiv = card.querySelector(".study-results");
  let awaitingRendered = false;

  const poll = async () => {
    let data;
    try {
      const res = await fetch(`/api/jobs/${jobId}`);
      data = await res.json();
    } catch (e) {
      stageText.textContent = "Lost contact with the server -- is it still running?";
      setTimeout(poll, 3000);
      return;
    }

    stageText.textContent = data.stage || "";
    logPanel.textContent = (data.log || []).join("\n");

    if (data.status === "running") {
      setTimeout(poll, 1500);
    } else if (data.status === "awaiting_active_space") {
      badge.textContent = "Confirm active space";
      badge.className = "badge running";
      if (!awaitingRendered) {
        awaitingRendered = true;
        renderStudyResults(resultsDiv, data.result);
        renderActiveSpaceConfirm(resultsDiv, jobId, data.result.active_space_suggestion, defaultCasscfNstate, false, () => {
          setTimeout(poll, 800);
        });
      }
      setTimeout(poll, 2000);
    } else if (data.status === "done") {
      stopBtn.remove();
      deleteBtn.disabled = false;
      badge.textContent = "Done";
      badge.className = "badge good";
      renderStudyResults(resultsDiv, data.result);
      // Only offered when this job actually went through CASSCF at least
      // once (i.e. an active-space suggestion exists) -- a plain RHF/CIS
      // job has nothing to add a combination to.
      if (data.result.active_space_suggestion) {
        renderTryAnotherButton(resultsDiv, jobId, data.result.active_space_suggestion, defaultCasscfNstate, () => {
          badge.textContent = "Running";
          badge.className = "badge running";
          deleteBtn.disabled = true;
          setTimeout(poll, 800);
        });
      }
    } else if (data.status === "cancelled") {
      stopBtn.remove();
      deleteBtn.disabled = false;
      badge.textContent = "Cancelled";
      badge.className = "badge cancelled";
    } else {
      stopBtn.remove();
      deleteBtn.disabled = false;
      badge.textContent = "Failed";
      badge.className = "badge bad";
      resultsDiv.innerHTML = `<p class="error">${data.error || "The job failed."}</p>`;
    }
  };
  poll();
}

function row(label, value) {
  return `<div class="result-row"><span class="label">${label}</span><span class="value">${value}</span></div>`;
}

// Per-stage rows show just the log file's own name -- the "Files" row
// above already gives the full out_dir, so repeating it on every log
// row is redundant clutter.
function basename(path) {
  return (path || "").split("/").pop();
}

function renderStudyResults(content, result) {
  if (!result) {
    content.innerHTML = "<p>No result data.</p>";
    return;
  }

  let html =
    row("Job", result.job_name) +
    row("RHF energy", `${result.rhf.energy.toFixed(6)} Hartree`) +
    row("NORB", result.rhf.norb) +
    row("RHF log", basename(result.rhf.log_path)) +
    row("Files", result.out_dir);

  if (result.cis) {
    let stateHead = "<th>State</th><th>Energy (Hartree)</th><th>S</th><th>Symmetry</th><th>Dominant transition</th>";
    let stateBody = "";
    result.cis.states.forEach((s) => {
      const top = s.transitions.reduce(
        (best, t) => (best === null || Math.abs(t.coefficient) > Math.abs(best.coefficient) ? t : best),
        null
      );
      const topDesc = top ? `${top.from_mo} &#8594; ${top.to_mo} (${top.coefficient.toFixed(4)})` : "n/a";
      stateBody += `<tr><td>${s.index}</td><td>${s.energy.toFixed(6)}</td><td>${s.spin}</td><td>${s.space_sym}</td><td>${topDesc}</td></tr>`;
    });
    html +=
      `<div class="results-table-wrap" style="margin-top:0.8rem"><table class="results-table"><thead><tr>${stateHead}</tr></thead><tbody>${stateBody}</tbody></table></div>` +
      row("CIS log", basename(result.cis.log_path));
  }

  content.innerHTML = html;

  // Each entry in result.combos is one independently-run active-space/state
  // combination (CASSCF -> optional XMCQDPT -> optional TRANSITN -> its own
  // energy table) -- per your request to compare several side by side, all
  // rendered here rather than just the most recent one.
  const combos = result.combos || [];
  if (combos.length > 1) {
    const allLatex = combos
      .filter((combo) => combo.energy_table)
      .map((combo) => combo.energy_table.latex)
      .join("\n\n");
    const combinedWrap = document.createElement("div");
    combinedWrap.style.marginTop = "0.8rem";
    combinedWrap.innerHTML = `<button type="button" class="small-btn combined-download-btn">Download all ${combos.length} tables (LaTeX)</button>`;
    content.appendChild(combinedWrap);
    combinedWrap.querySelector(".combined-download-btn").addEventListener("click", () => {
      downloadLatex(allLatex, "energy_tables.tex");
    });
  }

  combos.forEach((combo, i) => {
    const comboWrap = document.createElement("div");
    comboWrap.className = "describe-row";
    comboWrap.style.marginTop = "1rem";
    let comboHtml = `<label class="describe-label">Combination #${i + 1}</label>`;

    if (combo.error) {
      comboHtml += `<p class="error">${combo.error}</p>`;
      comboWrap.innerHTML = comboHtml;
      content.appendChild(comboWrap);
      return;
    }

    const c = combo.casscf;
    comboHtml += row("CASSCF", c.converged ? '<span class="badge good">converged</span>' : '<span class="badge bad">not converged</span>');
    if (c.converged) comboHtml += row("CASSCF final energy", `${c.final_energy.toFixed(6)} Hartree`);
    let stateHead = "<th>State</th><th>Energy (Hartree)</th>";
    let stateBody = "";
    c.state_energies.forEach(([index, energy]) => {
      stateBody += `<tr><td>${index}</td><td>${energy.toFixed(6)}</td></tr>`;
    });
    comboHtml +=
      `<div class="results-table-wrap" style="margin-top:0.8rem"><table class="results-table"><thead><tr>${stateHead}</tr></thead><tbody>${stateBody}</tbody></table></div>` +
      row("CASSCF log", basename(c.log_path)) +
      row("Active space", `NMCC=${c.active_space.nmcc} NDOC=${c.active_space.ndoc} NVAL=${c.active_space.nval}`);

    if (combo.xmcqdpt) {
      const x = combo.xmcqdpt;
      comboHtml += row("XMCQDPT", x.success ? '<span class="badge good">complete</span>' : '<span class="badge bad">did not succeed</span>');
      comboHtml += row("XMCQDPT log", basename(x.log_path));
    }

    if (combo.transitn) {
      const t = combo.transitn;
      comboHtml += row("Optical (TRANSITN, oscillator strengths)", t.success ? '<span class="badge good">complete</span>' : '<span class="badge bad">did not succeed</span>');
      comboHtml += row("Optical log", basename(t.log_path));
    }

    comboWrap.innerHTML = comboHtml;
    content.appendChild(comboWrap);

    if (combo.energy_table) {
      renderEnergyTable(comboWrap, combo.energy_table, `energy_table_${i + 1}.tex`);
    }
  });
}

function downloadLatex(latex, filename) {
  const blob = new Blob([latex], { type: "text/x-tex" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

function renderEnergyTable(content, table, downloadName) {
  const hasF = table.rows.some((r) => r.oscillator_strength != null);
  const wrap = document.createElement("div");
  wrap.style.marginTop = "0.8rem";
  let head = "<th>State</th><th>CASSCF (eV)</th><th>XMCQDPT (eV)</th>" + (hasF ? "<th>f</th>" : "");
  let body = "";
  table.rows.forEach((r) => {
    const xmcqdptStr = r.xmcqdpt_ev != null ? r.xmcqdpt_ev.toFixed(2) : "&ndash;";
    const fCell = hasF ? `<td>${r.oscillator_strength != null ? r.oscillator_strength.toFixed(4) : "&ndash;"}</td>` : "";
    body += `<tr><td>${r.label}</td><td>${r.casscf_ev.toFixed(2)}</td><td>${xmcqdptStr}</td>${fCell}</tr>`;
  });
  wrap.innerHTML = `
    <label class="describe-label">Vertical excitation energies (eV, relative to S0)</label>
    <div class="results-table-wrap"><table class="results-table"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>
    <div style="margin-top:0.6rem; display:flex; gap:0.5rem; align-items:center;">
      <button type="button" class="small-btn et-download-btn">Download LaTeX (.tex)</button>
      <button type="button" class="small-btn et-toggle-btn">View LaTeX source</button>
    </div>
    <textarea class="et-latex hidden" rows="12" readonly spellcheck="false" style="width:100%; margin-top:0.6rem; font-family:var(--mono); font-size:0.8rem; padding:0.65rem 0.75rem; background:var(--bg); color:var(--text); border:1px solid var(--panel-border); border-radius:6px;"></textarea>
  `;
  content.appendChild(wrap);

  const textarea = wrap.querySelector(".et-latex");
  textarea.value = table.latex;

  wrap.querySelector(".et-toggle-btn").addEventListener("click", () => {
    textarea.classList.toggle("hidden");
  });

  wrap.querySelector(".et-download-btn").addEventListener("click", () => {
    downloadLatex(table.latex, downloadName || "energy_table.tex");
  });
}

let _comboCardIdCounter = 0;

function _createComboCard(suggestion, defaultNstate, showMoSourceChoice) {
  const cardId = ++_comboCardIdCounter;
  const card = document.createElement("div");
  card.className = "as-combo-card";
  card.style.cssText = "border:1px solid var(--panel-border); border-radius:8px; padding:0.8rem; margin-top:0.7rem;";

  const orbitalCharacter = suggestion.orbital_character || {};
  const orbitalCharacterRows = [...suggestion.occ_selected, ...suggestion.virt_selected]
    .map((mo) => {
      const c = orbitalCharacter[String(mo)];
      return c
        ? `<tr><td>MO ${mo}</td><td>${c.label}</td><td>${c.coefficient.toFixed(4)}</td></tr>`
        : `<tr><td>MO ${mo}</td><td colspan="2">(not found in log)</td></tr>`;
    })
    .join("");
  const orbitalCharacterDetails = orbitalCharacterRows
    ? `
    <details class="as-orbital-character" style="margin-top:0.5rem">
      <summary>Orbital character (which atom/orbital each active MO is dominated by)</summary>
      <table class="results-table">
        <thead><tr><th>MO</th><th>Dominant atom/orbital</th><th>|coefficient|</th></tr></thead>
        <tbody>${orbitalCharacterRows}</tbody>
      </table>
    </details>`
    : "";

  // Only offered when there's actually a previous combination's
  // orbitals to draw from -- the very first combination in a job
  // always starts from the closed-shell (RHF) orbitals. Per your
  // request that reusing a previous combination's (possibly
  // MO-reordered) optimized orbitals be a choice, not automatic.
  const moSourceBlock = showMoSourceChoice
    ? `
    <div class="field" style="margin-top:0.5rem">
      <label>Starting orbitals for this combination's CASSCF</label>
      <div class="radio-group">
        <label><input type="radio" name="mo-source-${cardId}" class="as-mo-source" value="rhf" checked> Fresh closed-shell (RHF) orbitals</label>
        <label><input type="radio" name="mo-source-${cardId}" class="as-mo-source" value="previous"> Continue from the previous combination's optimized orbitals</label>
      </div>
    </div>`
    : "";

  card.innerHTML = `
    <div style="display:flex; justify-content:space-between; align-items:center;">
      <label class="describe-label" style="margin:0">Combination #<span class="as-card-index">1</span></label>
      <button type="button" class="small-btn as-remove-combo-btn" style="display:none">Remove</button>
    </div>
    ${orbitalCharacterDetails}
    <div class="field-row" style="margin-top:0.5rem">
      <div class="field">
        <label>Active occupied MOs (comma-separated)</label>
        <input type="text" class="as-occ" value="${suggestion.occ_selected.join(",")}">
      </div>
      <div class="field">
        <label>Active virtual MOs (comma-separated)</label>
        <input type="text" class="as-virt" value="${suggestion.virt_selected.join(",")}">
      </div>
    </div>
    <div class="field-row">
      <div class="field">
        <label>Number of CASSCF states (state-averaged)</label>
        <input type="number" class="as-nstate" value="${defaultNstate}" min="1">
      </div>
      <div class="field">
        <label>Memory (MWORDS) for CASSCF/XMCQDPT</label>
        <input type="number" class="as-mem-mwords" value="20" min="1">
      </div>
    </div>
    ${moSourceBlock}
    <label class="checkbox-row" style="margin-top:0.5rem"><input type="checkbox" class="as-xmcqdpt"> Also run XMCQDPT on top (dynamic correlation, same active space)</label>
    <label class="checkbox-row" style="margin-top:0.5rem"><input type="checkbox" class="as-transitn"> Also run oscillator strengths (TRANSITN, relative to S0)</label>
    <div class="as-recovery-section">
      <label class="checkbox-row" style="margin-top:0.5rem"><input type="checkbox" class="as-recovery"> If this doesn't converge, retry with more iterations, then fall back to a smaller active space and regrow</label>
      <div class="field-row as-recovery-fields hidden" style="margin-top:0.5rem">
        <div class="field">
          <label>Smaller active space: max electrons</label>
          <input type="number" class="as-smaller-electrons" value="${Math.min(4, 2 * suggestion.ndoc)}" min="2">
        </div>
        <div class="field">
          <label>Smaller active space: max orbitals</label>
          <input type="number" class="as-smaller-orbitals" value="${Math.min(4, suggestion.ndoc + suggestion.nval)}" min="2">
        </div>
      </div>
    </div>
  `;

  // The smaller-active-space recovery only makes sense when there's
  // room to shrink to something meaningfully smaller -- offering it for
  // an active space that's already (2,2) or similarly minimal is just
  // confusing clutter with nothing useful to fall back to, per your
  // report. Recomputed live off this card's own occ/virt fields.
  const recoverySection = card.querySelector(".as-recovery-section");
  const updateRecoveryVisibility = () => {
    const occCount = card.querySelector(".as-occ").value.split(",").map((x) => x.trim()).filter(Boolean).length;
    const virtCount = card.querySelector(".as-virt").value.split(",").map((x) => x.trim()).filter(Boolean).length;
    const applicable = (2 * occCount) >= 4 && (occCount + virtCount) >= 4;
    recoverySection.classList.toggle("hidden", !applicable);
  };
  updateRecoveryVisibility();
  card.querySelector(".as-occ").addEventListener("input", updateRecoveryVisibility);
  card.querySelector(".as-virt").addEventListener("input", updateRecoveryVisibility);
  card.querySelector(".as-recovery").addEventListener("change", () => {
    card.querySelector(".as-recovery-fields").classList.toggle("hidden", !card.querySelector(".as-recovery").checked);
  });

  return card;
}

function renderActiveSpaceConfirm(content, jobId, suggestion, defaultNstate, hasAnyPriorCombo, onSubmitted) {
  const box = document.createElement("div");
  box.className = "describe-row";
  const scoreLines = [...suggestion.occ_selected, ...suggestion.virt_selected]
    .map((mo) => `MO ${mo}: best |SAP coefficient| = ${(suggestion.scores[mo] || 0).toFixed(4)}`)
    .join("<br>");
  const cappedNote = suggestion.capped
    ? `<p class="hint" style="color:var(--bad)">Capped at the default max (16 electrons, 16 orbitals) -- dropped: occupied ${JSON.stringify(suggestion.occ_dropped)}, virtual ${JSON.stringify(suggestion.virt_dropped)}.</p>`
    : "";

  box.innerHTML = `
    <label class="describe-label">Suggested active space: NMCC=${suggestion.nmcc} NDOC=${suggestion.ndoc} NVAL=${suggestion.nval}
      (CAS(${2 * suggestion.ndoc} electrons, ${suggestion.ndoc + suggestion.nval} orbitals))</label>
    <p class="hint">${scoreLines}</p>
    ${cappedNote}
    <p class="hint">Add as many active-space/state combinations as you'd like to compare -- once you run them, they'll all execute one after another automatically (each through CASSCF, then XMCQDPT/TRANSITN if checked), with no further prompts.</p>
    <div class="as-combo-list"></div>
    <button type="button" class="small-btn as-add-combo-btn" style="margin-top:0.6rem">+ Add another combination</button>
    <div style="margin-top:0.9rem;">
      <button type="button" class="apply-btn as-run-all-btn">Run all combination(s)</button>
    </div>
    <p class="error as-error"></p>
  `;
  content.appendChild(box);

  const comboList = box.querySelector(".as-combo-list");

  const refreshCardChrome = () => {
    const cards = Array.from(comboList.querySelectorAll(".as-combo-card"));
    cards.forEach((card, i) => {
      card.querySelector(".as-card-index").textContent = i + 1;
      const removeBtn = card.querySelector(".as-remove-combo-btn");
      removeBtn.style.display = cards.length > 1 ? "" : "none";
    });
  };

  const addCombo = () => {
    // A card can offer "continue from the previous combination" whenever
    // there's actually a previous combination to draw from -- either one
    // already ran in an earlier batch (hasAnyPriorCombo), or an earlier
    // card in this same not-yet-submitted batch will have run by the
    // time this one executes.
    const showMoSourceChoice = hasAnyPriorCombo || comboList.children.length > 0;
    const card = _createComboCard(suggestion, defaultNstate, showMoSourceChoice);
    card.querySelector(".as-remove-combo-btn").addEventListener("click", () => {
      card.remove();
      refreshCardChrome();
    });
    comboList.appendChild(card);
    refreshCardChrome();
  };

  box.querySelector(".as-add-combo-btn").addEventListener("click", addCombo);
  addCombo();

  box.querySelector(".as-run-all-btn").addEventListener("click", async () => {
    const btn = box.querySelector(".as-run-all-btn");
    const errorEl = box.querySelector(".as-error");
    const cards = Array.from(comboList.querySelectorAll(".as-combo-card"));
    const combos = cards.map((card) => {
      const occ = card.querySelector(".as-occ").value.split(",").map((x) => parseInt(x.trim(), 10)).filter((x) => !isNaN(x));
      const virt = card.querySelector(".as-virt").value.split(",").map((x) => parseInt(x.trim(), 10)).filter((x) => !isNaN(x));
      const nstate = parseInt(card.querySelector(".as-nstate").value, 10) || 1;
      const memMwords = parseInt(card.querySelector(".as-mem-mwords").value, 10) || 20;
      const runXmcqdpt = card.querySelector(".as-xmcqdpt").checked;
      const runTransitn = card.querySelector(".as-transitn").checked;
      const recoveryApplicable = !card.querySelector(".as-recovery-section").classList.contains("hidden");
      const allowRecovery = recoveryApplicable && card.querySelector(".as-recovery").checked;
      const smallerElectrons = parseInt(card.querySelector(".as-smaller-electrons").value, 10) || 4;
      const smallerOrbitals = parseInt(card.querySelector(".as-smaller-orbitals").value, 10) || 4;
      const moSourceInput = card.querySelector(".as-mo-source:checked");
      const moSource = moSourceInput ? moSourceInput.value : "rhf";
      return {
        nstate, occ_selected: occ, virt_selected: virt, mem_mwords: memMwords, run_xmcqdpt: runXmcqdpt,
        run_transitn: runTransitn, allow_smaller_active_space_recovery: allowRecovery,
        smaller_max_electrons: smallerElectrons, smaller_max_orbitals: smallerOrbitals, mo_source: moSource,
      };
    });
    btn.disabled = true;
    errorEl.textContent = "";
    try {
      const res = await fetch(`/api/jobs/${jobId}/casscf`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ combos }),
      });
      if (!res.ok) {
        const body = await res.json().catch(() => ({}));
        throw new Error(body.detail || `HTTP ${res.status}`);
      }
      box.remove();
      onSubmitted();
    } catch (e) {
      btn.disabled = false;
      errorEl.textContent = `Couldn't start CASSCF: ${e.message}`;
    }
  });
}

function renderTryAnotherButton(content, jobId, suggestion, defaultNstate, onResumed) {
  const wrap = document.createElement("div");
  wrap.style.marginTop = "1rem";
  wrap.innerHTML = `<button type="button" class="small-btn try-another-btn">+ Try another combination</button>`;
  content.appendChild(wrap);
  wrap.querySelector(".try-another-btn").addEventListener("click", () => {
    wrap.remove();
    // Every combo card in this reopened batch can draw on a previous
    // combination's optimized orbitals -- at least one already ran, or
    // this button wouldn't be showing.
    renderActiveSpaceConfirm(content, jobId, suggestion, defaultNstate, true, onResumed);
  });
}
