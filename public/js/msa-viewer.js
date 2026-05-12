(function () {
  const MSA_SCRIPT_SRC = "https://cdn.jsdelivr.net/gh/thekaplanlab/msabrowser@master/javascript/msabrowser.js";

  let msaScriptPromise;

  function loadMsaScript() {
    if (window.MSABrowser && window.MSAProcessor) return Promise.resolve();
    if (msaScriptPromise) return msaScriptPromise;

    msaScriptPromise = new Promise((resolve, reject) => {
      const existing = document.querySelector(`script[src=\"${MSA_SCRIPT_SRC}\"]`);
      if (existing) {
        existing.addEventListener("load", () => resolve(), { once: true });
        existing.addEventListener("error", () => reject(new Error("Failed to load MSABrowser script.")), { once: true });
        return;
      }

      const script = document.createElement("script");
      script.src = MSA_SCRIPT_SRC;
      script.defer = true;
      script.onload = () => resolve();
      script.onerror = () => reject(new Error("Failed to load MSABrowser script."));
      document.head.appendChild(script);
    });

    return msaScriptPromise;
  }

  function parseJson(raw, fallback) {
    if (!raw) return fallback;
    try {
      return JSON.parse(raw);
    } catch (_error) {
      return fallback;
    }
  }

  function rewriteFastaHeadersForDisplay(fasta) {
    return fasta.replace(/^>([^\n]+)$/gm, (fullMatch, headerLine) => {
      const headerToken = headerLine.trim().split(/\s+/)[0];
      if (!headerToken.includes("|")) {
        return fullMatch;
      }

      return `>${headerToken.replace("|", " ")} ${headerLine.slice(headerToken.length).trim()}`.trimEnd();
    });
  }

  async function getFasta(el) {
    const fastaUrl = (el.dataset.fastaUrl || "").trim();
    const inlineFasta = (el.dataset.inlineFasta || "").trim();

    if (fastaUrl) {
      const response = await fetch(fastaUrl);
      if (!response.ok) throw new Error(`Failed to fetch FASTA from ${fastaUrl}`);
      return response.text();
    }

    if (inlineFasta) return inlineFasta;
    throw new Error("No FASTA data provided.");
  }

  async function initOne(el) {
    if (el.dataset.msaReady === "true") return;

    if (!el.id) {
      el.id = `msa-${Math.random().toString(36).slice(2, 10)}`;
    }

    try {
      await loadMsaScript();
      const fasta = await getFasta(el);
      const annotations = parseJson(el.dataset.annotations, []);
      const alterations = parseJson(el.dataset.alterations, []);
      const colorSchema = (el.dataset.colorSchema || "clustal").trim();
      const hasConsensus = el.dataset.hasConsensus !== "false";

      if (!window.__A2KOriginalMSAProcessor) {
        window.__A2KOriginalMSAProcessor = window.MSAProcessor;
        window.MSAProcessor = function wrappedMSAProcessor(options) {
          const originalFasta = options?.fasta || "";
          const viewerFasta = rewriteFastaHeadersForDisplay(originalFasta);
          const processed = window.__A2KOriginalMSAProcessor({
            ...options,
            fasta: viewerFasta,
          });
          processed.fasta = originalFasta;
          return processed;
        };
      }

      const viewer = new window.MSABrowser({
        id: el.id,
        msa: window.MSAProcessor({
          fasta,
          hasConsensus,
        }),
        annotations,
        alterations,
        colorSchema,
      });

      el.dataset.msaReady = "true";
      el._msaViewer = viewer;
    } catch (error) {
      el.dataset.msaReady = "error";
      // Keep errors visible in console for debugging bad data/paths.
      console.error("MSA viewer init failed", error);
    }
  }

  function initAll() {
    const targets = document.querySelectorAll("[data-msa]");
    targets.forEach((el) => {
      initOne(el);
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initAll, { once: true });
  } else {
    initAll();
  }

  document.addEventListener("astro:page-load", initAll);
})();
