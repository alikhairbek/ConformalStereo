# The program

## Start

- **Windows:** double-click `ConformalStereo` in the folder `ConformalStereo_App`.
- **macOS, Linux:** `python3 run.py` in this folder (or `sh start.sh`).

The program opens as a page in the web browser. It runs only on your computer: the page talks to the program that
you started, at the address 127.0.0.1, and to nothing else. Close the window of the program, or press *Quit* on the
page, to stop it; it also stops by itself a few minutes after the page was closed.

The program needs Python 3.8 or newer and no additional packages. On Windows the starter looks for a Python on the
computer. If there is none, it fetches the embeddable Python 3.12 from python.org into `ConformalStereo_App/runtime`
(about 11 MB, once). Delete the folder to remove everything. `ConformalStereo.bat --private-python` uses the private
copy even when another Python is installed.

## Pages

**Assign.** Give the DP4 or DP4+ result of one molecule and choose the cut-off.

| Input | What is read |
|---|---|
| DP4+ Excel sheet (`.xlsx`) | the stored values of the row "DP4+" for all isomers, with all digits, and the number of experimental shifts. Save the sheet in Excel after filling it in. |
| results of DP4+App (`.xlsx`) | the probabilities of the single nuclei in the sheets `p_sca` and `p_uns` (ten decimals), from which the unrounded likelihood is rebuilt. The sheet `results` keeps four decimals only. |
| output of PyDP4 or DP4-AI | the listed errors of the scaled shifts; the likelihood is computed again with the single Gaussian of PyDP4 (2.269 ppm for carbon, 0.187 ppm for proton), the model with which the DP4 workflows of the paper were calibrated. |
| typed values | one candidate per line, a name and a probability (per cent or fraction) or the logarithm of a likelihood. |

The number of digits sets the precision. If a value was rounded so far that a candidate can be neither kept nor
excluded (a printed 0.00%), the candidate is kept and marked.

The cut-off comes from one of the 76 workflows of the paper, from your own calibration, or from numbers you type.
A cut-off holds for the workflow with which it was calibrated: method, geometries, NMR level, solvent model and error
model. If your calculation differs in any of these, calibrate your own.

**Calibrate.** A table with one line per candidate of molecules of known structure:

    molecule,family,shifts,candidate,probability,actual
    compound 1,series A,24,isomer 1,97.3312,1
    compound 1,series A,24,isomer 2,2.6688,0

`family` is optional (molecules of one series or skeleton count together as one unit); `log_likelihood` may replace
`probability`. The 80% level needs at least 4 families, 90% needs 9 and 95% needs 19. A molecule that was read on the
page *Assign* can be added to the table with one click when its structure is known. The example on the page is one
workflow of the paper and gives its published cut-offs.

## Command line

    python conformalstereo.py assign result.xlsx --workflow "JOC2015 mPW1PW91/6-31+G**/PCM"
    python conformalstereo.py assign output.dp4 --calibration my_molecules.csv
    python conformalstereo.py assign --values "A=99.2134,B=0.7861,C=0.0005" --shifts 23 --cutoff 90=0.46
    python conformalstereo.py calibrate my_molecules.csv
    python conformalstereo.py workflows
    python conformalstereo.py selftest

`conformalstereo.py` can also be imported: `read_file`, `read_typed`, `read_table`, `calibrate`, `assign`, `cutoff`.

## Tests

    python tests/test_conformalstereo.py

28 tests, Python only. Inside the repository they also recompute the cut-offs of all 76 workflows from the pooled
table, compare the calibration with the tables of the analysis and repeat the worked example of Figure 1. The files in
`tests/fixtures` are made up; `tests/make_fixtures.py` writes them (the results file with the same calls as DP4+App 2.2).
Tested with Python 3.9, 3.11, 3.12 and 3.13.

## If something does not work

- *Windows protected your PC*: the file came from the internet. Choose *More info*, then *Run anyway*.
- *Python could not be fetched*: the computer could not reach python.org. Install Python from python.org and start again.
- The page does not open: copy the address that the window of the program shows into the browser.
- A workbook is refused: files in the old `.xls` format have to be saved as `.xlsx`; the DP4+ sheet has to be saved
  in Excel after the data were entered, so that it stores its results.
