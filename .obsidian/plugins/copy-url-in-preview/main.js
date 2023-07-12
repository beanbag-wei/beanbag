/*
THIS IS A GENERATED/BUNDLED FILE BY ESBUILD
if you want to view the source, please visit the github repository of this plugin
*/

var __defProp = Object.defineProperty;
var __getOwnPropDesc = Object.getOwnPropertyDescriptor;
var __getOwnPropNames = Object.getOwnPropertyNames;
var __hasOwnProp = Object.prototype.hasOwnProperty;
var __export = (target, all) => {
  for (var name in all)
    __defProp(target, name, { get: all[name], enumerable: true });
};
var __copyProps = (to, from, except, desc) => {
  if (from && typeof from === "object" || typeof from === "function") {
    for (let key of __getOwnPropNames(from))
      if (!__hasOwnProp.call(to, key) && key !== except)
        __defProp(to, key, { get: () => from[key], enumerable: !(desc = __getOwnPropDesc(from, key)) || desc.enumerable });
  }
  return to;
};
var __toCommonJS = (mod) => __copyProps(__defProp({}, "__esModule", { value: true }), mod);

// src/main.ts
var main_exports = {};
__export(main_exports, {
  default: () => CopyUrlInPreview
});
module.exports = __toCommonJS(main_exports);
var import_obsidian = require("obsidian");

// src/helpers.ts
var loadImageBlobTimeout = 5e3;
function withTimeout(ms, promise) {
  const timeout = new Promise((resolve, reject) => {
    const id = setTimeout(() => {
      clearTimeout(id);
      reject(`timed out after ${ms} ms`);
    }, ms);
  });
  return Promise.race([
    promise,
    timeout
  ]);
}
async function loadImageBlob(imgSrc) {
  const loadImageBlobCore = () => {
    return new Promise((resolve, reject) => {
      const image = new Image();
      image.crossOrigin = "anonymous";
      image.onload = () => {
        const canvas = document.createElement("canvas");
        canvas.width = image.width;
        canvas.height = image.height;
        const ctx = canvas.getContext("2d");
        ctx.drawImage(image, 0, 0);
        canvas.toBlob((blob) => {
          resolve(blob);
        });
      };
      image.onerror = async () => {
        try {
          await fetch(image.src, { "mode": "no-cors" });
          const blob = await loadImageBlob(`https://api.allorigins.win/raw?url=${encodeURIComponent(imgSrc)}`);
          resolve(blob);
        } catch (e) {
          reject();
        }
      };
      image.src = imgSrc;
    });
  };
  return withTimeout(loadImageBlobTimeout, loadImageBlobCore());
}
function onElement(el, event, selector, listener, options) {
  el.on(event, selector, listener, options);
  return () => el.off(event, selector, listener, options);
}

// src/main.ts
var IMAGE_URL_PREFIX = "/_capacitor_file_";
var SUCCESS_NOTICE_TIMEOUT = 1800;
var longTapTimeout = 500;
var deleteTempFileTimeout = 6e4;
var OPEN_PDF_MENU_BORDER_SIZE = 100;
var OPEN_PDF_MENU_TIMEOUT = 5e3;
var CopyUrlInPreview = class extends import_obsidian.Plugin {
  constructor() {
    super(...arguments);
    this.longTapTimeoutId = null;
  }
  onload() {
    this.registerDocument(document);
    app.workspace.on(
      "window-open",
      (workspaceWindow, window2) => {
        this.registerDocument(window2.document);
      }
    );
  }
  registerDocument(document2) {
    this.register(
      onElement(
        document2,
        "contextmenu",
        "a.external-link",
        this.onClick.bind(this)
      )
    );
    this.register(
      onElement(
        document2,
        "mouseover",
        ".pdf-embed iframe, .pdf-embed div.pdf-container, .workspace-leaf-content[data-type=pdf]",
        this.showOpenPdfMenu.bind(this)
      )
    );
    this.register(
      onElement(
        document2,
        "mousemove",
        ".pdf-canvas",
        this.showOpenPdfMenu.bind(this)
      )
    );
    if (import_obsidian.Platform.isDesktop) {
      this.register(
        onElement(
          document2,
          "contextmenu",
          "img",
          this.onClick.bind(this)
        )
      );
      this.register(
        onElement(
          document2,
          "mouseover",
          ".cm-link, .cm-hmd-internal-link",
          this.storeLastHoveredLinkInEditor.bind(this)
        )
      );
      this.register(
        onElement(
          document2,
          "mouseover",
          "a.internal-link",
          this.storeLastHoveredLinkInPreview.bind(this)
        )
      );
    } else {
      this.register(
        onElement(
          document2,
          "touchstart",
          "img",
          this.startWaitingForLongTap.bind(this)
        )
      );
      this.register(
        onElement(
          document2,
          "touchend",
          "img",
          this.stopWaitingForLongTap.bind(this)
        )
      );
      this.register(
        onElement(
          document2,
          "touchmove",
          "img",
          this.stopWaitingForLongTap.bind(this)
        )
      );
    }
  }
  storeLastHoveredLinkInEditor(event) {
    var _a;
    const editor = (_a = app.workspace.getActiveViewOfType(import_obsidian.MarkdownView)) == null ? void 0 : _a.editor;
    if (!editor) {
      return;
    }
    const position = editor.posAtMouse(event);
    const token = editor.getClickableTokenAt(position);
    if (!token) {
      return;
    }
    this.lastHoveredLinkTarget = token.text;
  }
  storeLastHoveredLinkInPreview(event, link) {
    this.lastHoveredLinkTarget = link.getAttribute("data-href");
  }
  showOpenPdfMenu(event, el) {
    var _a;
    if (this.openPdfMenu || this.preventReopenPdfMenu) {
      return;
    }
    const rect = el.getBoundingClientRect();
    if (rect.left + OPEN_PDF_MENU_BORDER_SIZE < event.x && event.x < rect.right - OPEN_PDF_MENU_BORDER_SIZE && rect.top + OPEN_PDF_MENU_BORDER_SIZE < event.y && event.y < rect.bottom - OPEN_PDF_MENU_BORDER_SIZE) {
      return;
    }
    const pdfEmbed = el.closest(".pdf-embed");
    let pdfFile;
    if (pdfEmbed) {
      let pdfLink;
      if (pdfEmbed.hasClass("popover")) {
        pdfLink = this.lastHoveredLinkTarget;
      } else {
        pdfLink = (_a = pdfEmbed.getAttr("src")) != null ? _a : this.lastHoveredLinkTarget;
      }
      pdfLink = pdfLink == null ? void 0 : pdfLink.replace(/#page=\d+$/, "");
      const currentNotePath = this.app.workspace.getActiveFile().path;
      pdfFile = this.app.metadataCache.getFirstLinkpathDest(pdfLink, currentNotePath);
    } else {
      pdfFile = this.app.workspace.getActiveFile();
    }
    const menu = new import_obsidian.Menu();
    this.registerEscapeButton(menu);
    menu.onHide(() => this.openPdfMenu = null);
    menu.addItem(
      (item) => item.setIcon("pdf-file").setTitle("Open PDF externally").onClick(async () => {
        this.preventReopenPdfMenu = true;
        setTimeout(() => {
          this.preventReopenPdfMenu = false;
        }, OPEN_PDF_MENU_TIMEOUT);
        this.hideOpenPdfMenu();
        if (import_obsidian.Platform.isDesktop) {
          await this.app.openWithDefaultApp(pdfFile.path);
        } else {
          await this.app.vault.adapter.open(pdfFile.path);
        }
      })
    );
    menu.showAtMouseEvent(event);
    this.openPdfMenu = menu;
    setTimeout(this.hideOpenPdfMenu.bind(this), OPEN_PDF_MENU_TIMEOUT);
  }
  registerEscapeButton(menu, document2 = activeDocument) {
    menu.register(
      onElement(
        document2,
        "keydown",
        "*",
        (e) => {
          if (e.key === "Escape") {
            e.preventDefault();
            e.stopPropagation();
            menu.hide();
          }
        }
      )
    );
  }
  hideOpenPdfMenu() {
    if (this.openPdfMenu) {
      this.openPdfMenu.hide();
    }
  }
  // mobile
  startWaitingForLongTap(event, img) {
    if (this.longTapTimeoutId) {
      clearTimeout(this.longTapTimeoutId);
      this.longTapTimeoutId = null;
    } else {
      if (event.targetTouches.length == 1) {
        this.longTapTimeoutId = window.setTimeout(this.processLongTap.bind(this, event, img), longTapTimeout);
      }
    }
  }
  // mobile
  stopWaitingForLongTap() {
    if (this.longTapTimeoutId) {
      clearTimeout(this.longTapTimeoutId);
      this.longTapTimeoutId = null;
    }
  }
  // mobile
  async processLongTap(event, img) {
    event.stopPropagation();
    this.longTapTimeoutId = null;
    const adapter = this.app.vault.adapter;
    const electronWindow = window;
    const basePath = adapter.getFullPath("");
    const webviewServerUrl = electronWindow.WEBVIEW_SERVER_URL;
    const localImagePrefixUrl = webviewServerUrl + IMAGE_URL_PREFIX + basePath;
    if (img.src.startsWith(localImagePrefixUrl)) {
      const encodedImageFileRelativePath = img.src.replace(localImagePrefixUrl, "");
      const imageFileRelativePath = decodeURIComponent(encodedImageFileRelativePath);
      await adapter.open(imageFileRelativePath);
    } else {
      try {
        const blob = await loadImageBlob(img.src);
        if (!blob.type.startsWith("image/")) {
          new import_obsidian.Notice(`Unsupported mime type ${blob.type}`);
          return;
        }
        const extension = blob.type.replace("image/", "");
        const randomGuid = window.URL.createObjectURL(new Blob([])).split("/").pop();
        const tempFileName = `/.temp-${randomGuid}.${extension}`;
        const buffer = await blob.arrayBuffer();
        await adapter.writeBinary(tempFileName, buffer);
        setTimeout(() => adapter.remove(tempFileName), deleteTempFileTimeout);
        new import_obsidian.Notice("Image was temporarily saved and will be removed in 1 minute");
        await adapter.open(tempFileName);
      } catch (e) {
        new import_obsidian.Notice("Cannot open image");
      }
    }
  }
  // Android gives a PointerEvent, a child to MouseEvent.
  // Positions are not accurate from PointerEvent.
  // There's also TouchEvent
  // The event has target, path, toEvent (null on Android) for finding the link
  onClick(event) {
    event.preventDefault();
    const target = event.target;
    const imgType = target.localName;
    const menu = new import_obsidian.Menu();
    switch (imgType) {
      case "img": {
        const image = target.currentSrc;
        const url = new URL(image);
        const protocol = url.protocol;
        switch (protocol) {
          case "app:":
          case "data:":
          case "http:":
          case "https:":
            menu.addItem(
              (item) => item.setIcon("image-file").setTitle("Copy image to clipboard").onClick(async () => {
                try {
                  const blob = await loadImageBlob(image);
                  const data = new ClipboardItem({ [blob.type]: blob });
                  await navigator.clipboard.write([data]);
                  new import_obsidian.Notice("Image copied to the clipboard!", SUCCESS_NOTICE_TIMEOUT);
                } catch (e) {
                  new import_obsidian.Notice("Error, could not copy the image!");
                }
              })
            );
            if (protocol === "app:" && import_obsidian.Platform.isDesktop) {
              const baseFilePath = app.vault.adapter.getFilePath("");
              const baseFilePathName = baseFilePath.pathname;
              const urlPathName = url.pathname;
              if (urlPathName.startsWith(baseFilePathName)) {
                let relativePath = urlPathName.replace(baseFilePathName, "");
                relativePath = decodeURI(relativePath);
                menu.addItem(
                  (item) => item.setIcon("arrow-up-right").setTitle("Open in default app").onClick(() => app.openWithDefaultApp(relativePath))
                );
                menu.addItem(
                  (item) => item.setIcon("arrow-up-right").setTitle(import_obsidian.Platform.isMacOS ? "Reveal in finder" : "Show in system explorer").onClick(() => {
                    app.showInFolder(relativePath);
                  })
                );
                menu.addItem(
                  (item) => item.setIcon("folder").setTitle("Reveal file in navigation").onClick(() => {
                    const abstractFilePath = app.vault.getAbstractFileByPath(relativePath.substring(1));
                    app.internalPlugins.getEnabledPluginById("file-explorer").revealInFolder(abstractFilePath);
                  })
                );
              }
            }
            break;
          default:
            new import_obsidian.Notice(`no handler for ${protocol} protocol`);
            return;
        }
        break;
      }
      case "a": {
        const link = target.href;
        menu.addItem(
          (item) => item.setIcon("link").setTitle("Copy URL").onClick(() => {
            navigator.clipboard.writeText(link);
            new import_obsidian.Notice("URL copied to your clipboard", SUCCESS_NOTICE_TIMEOUT);
          })
        );
        break;
      }
      default:
        new import_obsidian.Notice("No handler for this image type!");
        return;
    }
    this.registerEscapeButton(menu);
    menu.showAtPosition({ x: event.pageX, y: event.pageY });
    this.app.workspace.trigger("copy-url-in-preview:contextmenu", menu);
  }
};