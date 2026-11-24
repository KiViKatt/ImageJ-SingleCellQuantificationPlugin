# Landing page

One static page that introduces Cell Quant Studio to people outside the field. There is nothing to build: open `index.html` in a browser.

| File or folder | What it is |
|---|---|
| `index.html` | The whole page: text, styles and scripts in one file. |
| `img/` | Pictures of the application window, recorded from the real window (see below). |
| `fonts/` | Chakra Petch and IBM Plex Sans with their licence texts (SIL Open Font License 1.1). |
| `record/` | The scripts that made the pictures. |
| `_headers` | For Cloudflare: two safety headers for every file and a long cache for the typefaces. |

## What is on the page

- **A light and a dark theme.** The page follows the system; the round button in the top bar switches, and the choice is remembered in the browser. The colors are those of the application's Light and Dark looks.
- **The demo.** Seven steps, following the one-click demo in [`../demo/`](../demo/): open the image, find the tissue, analyse (the Researcher view, which the Results tab shows first), see the map (the Pitch view), point at a zone, switch the marker, check the picture. It starts by itself when it is scrolled into view and has a Play and Pause button; picking a step stops it. Without scripts the step buttons still work.
- **The ring map at the top** shows the numbers of the demo run.
- **A QR code** for `https://quantmycells.com`, beside the two buttons on the first screen and again in the install panel, for people who watch the page on a shared screen or a stream. It is always black on white, and it is hidden on phones. The card on the first screen is exactly as tall as the two buttons.
- **Before and after.** Plain Fiji with its own tools on the left, Cell Quant Studio on the right, for starting an analysis and for reading the answer; under it, the chores the application takes over.

## The pictures are from a synthetic image

The demo was recorded on a test image made by `record/make_demo_image.py`, not on real tissue. Nuclei are drawn as blurred ovals, and a share of them, larger near the surface, is lit in each marker channel, so that the depth zones differ. The page says so under the demo. Its numbers say nothing about biology.

In that run 2,600 nuclei were drawn and the headline count is 2,699. The difference is expected: the three depth zones share borders, a nucleus lying on a border is clipped there, and each part that passes the size filter is counted in its own zone.

## The "before" pictures

`img/before-1.jpg` and `img/before-2.jpg` show plain Fiji doing the same job on the same synthetic image. `record/plain_fiji.py` opens Fiji's own windows for it (the image with a threshold, the Channels Tool, the Threshold window, ROI Manager, the Analyze Particles dialog; then the Results and Summary tables). Each window was captured by itself and the captures were arranged on a plain background, so the pictures show real Fiji windows but not one real screen. The page says so under them.

## The QR code

The code is drawn in `index.html` as an inline SVG (the two `svg class="qr"` elements): 25 x 25 modules, error correction level M, with a quiet margin of four modules. If the address changes, make a new code for the new address, replace both drawings, and check with a phone that it opens the right page.

## Recording the pictures again

Do this after the window changes. From the repository folder, with the Fiji folder in place of `<Fiji>`:

```text
mkdir -p /data/organoids

CQ_DEMO_IMAGE="/data/organoids/Organoid section (synthetic).ome.tif" \
CELLQUANT_JAVA_OPTS=-Xmx8g \
bash tests/run_tests.sh <Fiji> "$PWD/site/record/make_demo_image.py"

CQ_DEMO_IMAGE="/data/organoids/Organoid section (synthetic).ome.tif" \
CQ_RESULTS=/data/organoids/results \
CQ_SHOTS=/data/organoids/shots \
CELLQUANT_JAVA_OPTS=-Xmx10g \
bash tests/run_tests.sh <Fiji> "$PWD/site/record/record_demo.py"
```

- Use a folder whose path you do not mind publishing: the image path and the output folder are shown in the window, so they are in the pictures.
- Do not click in the application window while it records.
- The pictures arrive as PNG files, 2880 x 1744. The page uses them as JPEG, 2160 pixels wide (the zones picture 1308 wide):

| Recorded | On the page |
|---|---|
| `dark_2_profile.png`, `light_2_profile.png` | `img/demo-1.jpg`, `img/demo-1-light.jpg` |
| `zones.png` | `img/demo-2.jpg` |
| `dark_3_research.png`, `light_3_research.png` | `img/demo-3.jpg`, `img/demo-3-light.jpg` |
| `dark_3_results.png`, `light_3_results.png` | `img/demo-4.jpg`, `img/demo-4-light.jpg` |
| `dark_4_point.png`, `light_4_point.png` | `img/demo-5.jpg`, `img/demo-5-light.jpg` |
| `dark_5_marker.png`, `light_5_marker.png` | `img/demo-6.jpg`, `img/demo-6-light.jpg` |
| `dark_6_picture.png`, `light_6_picture.png` | `img/demo-7.jpg`, `img/demo-7-light.jpg` |

If the numbers of the run change, change them in the ring map at the top of `index.html` too (the total, the three zone values and the lengths of the arcs). The pictures now on the page were recorded before the noise of the test image was seeded, so the next recording can differ from them by a few nuclei; after that the image, and with it the numbers, stay the same.

## Publishing

The page is published on Cloudflare with its `wrangler` tool, as plain static files. Copy the page's files, without this README and `record/`, to a clean folder and deploy that:

```text
rsync -a --exclude README.md --exclude record --exclude .DS_Store site/ /tmp/site_dist/
npx wrangler deploy /tmp/site_dist --name quantmycells --compatibility-date 2026-11-19 \
    --domain quantmycells.com --domain www.quantmycells.com
```

`npx wrangler login` signs in first; it opens Cloudflare in the browser. The two `--domain` options attach the page to its address (the domain must be in the same Cloudflare account); with them the `workers.dev` address is switched off. Run the command again after any change to publish the new version.

Any other static host works too: point it at the same files, with no build command. Nothing on the page calls another server.

## What was checked

- Layout without overflow from 320 to 1920 pixels wide, in both themes (headless Chrome).
- The demo player, the before and after switch, the fade-in of the sections and the count in the ring map, in Chromium at a laptop and a phone size.
- The QR code: read back from pictures of the first screen at window sizes from 1000 x 700 to 2560 x 1300, and from the same pictures shrunk and compressed like a 720p stream. Not tried with a phone camera on a real stream.
- On the published address: every file the page asks for loads, and a wrong path gives a 404.
- Safari: the layout was looked at once, at a laptop and a narrow width, before the QR code and the before and after section were added. Nothing was clicked there.
- Not checked: Firefox, and a real phone.
