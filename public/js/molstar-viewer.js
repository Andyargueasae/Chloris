(function () {
  let styleInjected = false;

  function injectSpinStyles() {
    if (styleInjected) return;
    styleInjected = true;

    const style = document.createElement("style");
    style.textContent = "@keyframes a2k-spin-viewer { from { transform: rotate(0deg); } to { transform: rotate(360deg); } }" +
      "[data-molstar].is-spinning canvas { transform-origin: center center; animation: a2k-spin-viewer 8s linear infinite; }";
    document.head.appendChild(style);
  }

  function detectFormat(path) {
    const lower = (path || "").toLowerCase();
    if (lower.endsWith(".bcif") || lower.endsWith(".bcif.gz")) return { format: "bcif", isBinary: true };
    if (lower.endsWith(".pdb")) return { format: "pdb", isBinary: false };
    return { format: "mmcif", isBinary: false };
  }

  async function waitForMolstar() {
    for (let i = 0; i < 120; i += 1) {
      if (window.molstar && window.molstar.Viewer && typeof window.molstar.Viewer.create === "function") {
        return window.molstar;
      }
      await new Promise((resolve) => setTimeout(resolve, 50));
    }
    throw new Error("Mol* global not available. Ensure /vendor/molstar/molstar.js loads before /js/molstar-viewer.js.");
  }

  async function initOne(container) {
    if (container.dataset.molstarReady === "true") return;

    const section = container.closest("[data-protein-viewer]");
    if (!section) return;

    try {
      injectSpinStyles();
      const molstar = await waitForMolstar();

      const viewer = await molstar.Viewer.create(container, {
        layoutShowControls: true,
        layoutShowSequence: true,
        layoutShowLog: false,
        layoutIsExpanded: false,
      });

      const cifUrl = (container.dataset.cifUrl || "").trim();
      const background = (container.dataset.background || "dark").trim();
      const shouldReset = container.dataset.resetCamera !== "false";

      if (background === "light") {
        container.style.background = "#ffffff";
      }

      if (cifUrl) {
        const guessed = detectFormat(cifUrl);
        await viewer.loadStructureFromUrl(cifUrl, guessed.format, guessed.isBinary);
      }

      const resetBtn = section.querySelector(".js-reset");
      const screenshotBtn = section.querySelector(".js-screenshot");
      const spinBtn = section.querySelector(".js-spin");
      const fullscreenBtn = section.querySelector(".js-fullscreen");
      const fileInput = section.querySelector(".js-file");

      if (resetBtn) {
        resetBtn.addEventListener("click", () => {
          viewer.plugin.canvas3d?.requestCameraReset?.({ durationMs: 300 });
        });
      }

      if (screenshotBtn) {
        screenshotBtn.addEventListener("click", () => {
          const helper = viewer.plugin.helpers?.viewportScreenshot;
          if (helper?.download) {
            helper.download();
            return;
          }

          const canvas = container.querySelector("canvas");
          if (!canvas) return;
          const a = document.createElement("a");
          a.href = canvas.toDataURL("image/png");
          a.download = "structure.png";
          a.click();
        });
      }

      if (spinBtn) {
        spinBtn.addEventListener("click", () => {
          container.classList.toggle("is-spinning");
        });
      }

      if (fullscreenBtn) {
        fullscreenBtn.addEventListener("click", () => {
          section.classList.toggle("is-fullscreen");
          fullscreenBtn.classList.toggle("is-active");
          
          // Close fullscreen with Escape key
          if (section.classList.contains("is-fullscreen")) {
            const handleEscape = (e) => {
              if (e.key === "Escape") {
                section.classList.remove("is-fullscreen");
                fullscreenBtn.classList.remove("is-active");
                document.removeEventListener("keydown", handleEscape);
              }
            };
            document.addEventListener("keydown", handleEscape);
          }
        });
      }

      if (fileInput) {
        fileInput.addEventListener("change", async (event) => {
          const file = event.target.files && event.target.files[0];
          if (!file) return;

          try {
            const guessed = detectFormat(file.name);
            if (guessed.isBinary) {
              const data = new Uint8Array(await file.arrayBuffer());
              await viewer.loadStructureFromData(data, guessed.format, { dataLabel: file.name });
            } else {
              const text = await file.text();
              await viewer.loadStructureFromData(text, guessed.format, { dataLabel: file.name });
            }

            if (shouldReset) {
              viewer.plugin.canvas3d?.requestCameraReset?.({ durationMs: 200 });
            }
          } catch (error) {
            console.error("Failed to load local structure", error);
          }
        });
      }

      if (shouldReset) {
        viewer.plugin.canvas3d?.requestCameraReset?.({ durationMs: 200 });
      }

      container.dataset.molstarReady = "true";
      container._molstarViewer = viewer;
    } catch (error) {
      container.dataset.molstarReady = "error";
      console.error("Mol* viewer init failed", error);
    }
  }

  function initAll() {
    const targets = document.querySelectorAll("[data-molstar]");
    targets.forEach((container) => {
      initOne(container);
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", initAll, { once: true });
  } else {
    initAll();
  }

  document.addEventListener("astro:page-load", initAll);
})();
