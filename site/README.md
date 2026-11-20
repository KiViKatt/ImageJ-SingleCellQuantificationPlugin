# Landing page

One static page that introduces Cell Quant Studio to people outside the field. There is nothing to build: open `index.html` in a browser.

| File or folder | What it is |
|---|---|
| `index.html` | The whole page: text, styles and scripts in one file. |
| `img/` | Pictures of the application window, recorded from the real window (see below). |
| `fonts/` | Chakra Petch and IBM Plex Sans with their licence texts (SIL Open Font License 1.1). |
| `record/` | The two scripts that made the pictures. |

## What is on the page

- **A light and a dark theme.** The page follows the system; the round button in the top bar switches, and the choice is remembered in the browser. The colors are those of the application's Light and Dark looks.
- **The demo.** Six steps, the same as the one-click demo in [`../demo/`](../demo/): open the image, find the tissue, analyse, point at a zone, switch the marker, check the picture. It starts by itself when it is scrolled into view and has a Play and Pause button; picking a step stops it. Without scripts the step buttons still work.
- **The ring map at the top** shows the numbers of the demo run.

## The pictures are from a synthetic image

The demo was recorded on a test image made by `record/make_demo_image.py`, not on real tissue. Nuclei are drawn as blurred ovals, and a share of them, larger near the surface, is lit in each marker channel, so that the depth zones differ. The page says so under the demo. Its numbers say nothing about biology.

In that run 2,600 nuclei were drawn and the headline count is 2,703. The difference is expected: the three depth zones share borders, a nucleus lying on a border is clipped there, and each part that passes the size filter is counted in its own zone.

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
| `dark_3_results.png`, `light_3_results.png` | `img/demo-3.jpg`, `img/demo-3-light.jpg` |
| `dark_4_point.png`, `light_4_point.png` | `img/demo-4.jpg`, `img/demo-4-light.jpg` |
| `dark_5_marker.png`, `light_5_marker.png` | `img/demo-5.jpg`, `img/demo-5-light.jpg` |
| `dark_6_picture.png`, `light_6_picture.png` | `img/demo-6.jpg`, `img/demo-6-light.jpg` |

If the numbers of the run change, change them in the ring map at the top of `index.html` too (the total, the three zone values and the lengths of the arcs).

## Publishing

Any static host works. Point it at this folder (`site`) with no build command. Nothing on the page calls another server.

## What was checked

- Layout without overflow or overlap from 320 to 1920 pixels wide, in both themes (headless Chrome).
- The demo player, the fade-in of the sections and the count in the ring map, in Chromium at a laptop and a phone size.
- Not checked: Safari, Firefox, and a real phone.
