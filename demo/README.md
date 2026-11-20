# One-click demo

`Cell_Quant_Studio_Demo.py` shows Cell Quant Studio end to end, in about one minute. It:

1. opens the application in the Dark look;
2. loads a demo image and a demo settings profile;
3. finds the tissue outline;
4. adds three depth zones that do not overlap: 0-100 µm, 100-200 µm, and over 200 µm, measured from the tissue edge, not from the edge of the picture;
5. runs the analysis;
6. ends on the Results tab.

It changes nothing in how images are measured. It only clicks through the window the way a user would.

## Set up once

1. Copy `Cell_Quant_Studio_Demo.py` beside `Cell_Quant_Studio_final.py` in `Fiji.app/plugins/Scripts/Cell Quant Studio/`. It then appears in Fiji's Plugins menu.
2. Copy a profile there too. `examples/Demo_4Channel_HOPX_PAX6_DAPI_TBR2.properties` fits files with HOPX in C1, PAX6 in C2, DAPI in C3 and TBR2 in C4.
3. Optional: a file `Cell_Quant_Studio_Demo.txt` beside the script, to start without any questions:

   ```text
   image=C:/path/to/demo image.czi
   profile=C:/path/to/Fiji.app/plugins/Scripts/Cell Quant Studio/Demo_profile.properties
   look=Dark
   ```

   Without it, the first run asks for the image and the profile and remembers them. Hold Shift while starting the demo to choose again.

Results go to `C:\CQ_demo` (on macOS and Linux, `~/CQ_demo`).

## Run sheet for a 3-minute presentation

| When | Show | Say |
|---|---|---|
| Before | Fiji open, nothing else | "Counting cells in organoid sections is usually done by eye, threshold by threshold, and nobody can rerun it." |
| 0:00 | Plugins > Cell Quant Studio > Cell Quant Studio Demo | "One click: the same steps a scientist takes, recorded so anyone can repeat them." |
| 0:10 | Session tab: the image, its size in µm, the channel list | "It reads the microscope file: every channel with its wavelength. A channel called DAPI that isn't the 405 nm laser is flagged. That mistake happens in real labs." |
| 0:25 | Tissue outline and the three zones | "It finds the organoid outline and divides it by depth from the surface." |
| 0:35 | The analysis runs | "Thousands of nuclei counted, four markers measured, every setting saved next to the numbers." |
| 0:50 | Results: the rings | "Point at a ring: each zone lights up. Outer zone, middle, core." |
| 1:10 | Click "More choices", choose TBR2 or HOPX | "The neuronal progenitor markers sit near the surface and fade toward the core: the layering of a developing cortex." |
| 1:40 | Click the picture | "Every counted nucleus is outlined. You can check it, not just trust it." |
| 2:00 | Light / Dark switch | "Light or dark, the same tool." |
| 2:15 | "Show result files" | "Plain tables any lab can open, plus the exact settings. Open source, built in Fiji, by scientists for scientists." |

**If something goes wrong:** the Results tab keeps the last finished analysis. Click the Results tab and present from there.
