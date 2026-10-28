# GTP Intersection Numbers

Companion code for the paper

> H. Argüz, A. Bourget, P. Bousseau, J. Grimminger, S. Schäfer-Nameki,
> *5d SCFTs and Intersection Theory for Generalized Toric Polygons* (2026).

Given a generalized toric polygon (GTP) with a unimodular triangulation, the code computes the triple intersection numbers of the divisors of the associated Calabi–Yau threefold. 

## Requirements

[Wolfram Mathematica](https://www.wolfram.com/mathematica/) with its notebook interface (tested with version 14.3). The interactive editor needs the notebook front end and does not run from `wolframscript`.

## Usage

1. Download `GTP_Intersection_Numbers.nb` and `GTPIntersections.wl` into the same folder. The notebook loads the package from its own folder.
2. Open `GTP_Intersection_Numbers.nb` and evaluate the section **Load the package**. If Mathematica asks whether to evaluate the initialization cells, answer *Yes*.
3. Go to the section **Draw a GTP and compute everything** and evaluate `res = gtpRun[];`.
   - Enter the size of the grid.
   - In the window that opens, draw the edges of the triangulation and the 7-brane cuts.
   - Click **Done**.

   The output is the GTP, its triple intersection numbers and the TikZ code of the picture.

The section **Examples of Section 7** reproduces the examples of the paper.

## Citation

If you use this code, please cite the paper above.
