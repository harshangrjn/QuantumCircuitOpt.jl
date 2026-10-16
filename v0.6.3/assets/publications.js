document.addEventListener("DOMContentLoaded", () => {
  const copyText = async (text) => {
    if (navigator.clipboard && window.isSecureContext) {
      await navigator.clipboard.writeText(text);
      return;
    }

    const textarea = document.createElement("textarea");
    textarea.value = text;
    textarea.style.position = "fixed";
    textarea.style.opacity = "0";
    document.body.appendChild(textarea);
    textarea.select();
    document.execCommand("copy");
    textarea.remove();
  };

  document.querySelectorAll(".bibtex-copy").forEach((button) => {
    button.addEventListener("click", async () => {
      const citation = button
        .closest(".citation-block")
        .querySelector("code").textContent;

      try {
        await copyText(citation);
        button.textContent = "Copied";
      } catch {
        button.textContent = "Copy failed";
      }

      window.setTimeout(() => {
        button.textContent = "Copy BibTeX";
      }, 1600);
    });
  });
});
